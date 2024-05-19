import asyncio
from collections import defaultdict

import pandas as pd

from metamist.apis import AssayApi
from metamist.graphql import gql, query
from metamist.models import AssayUpsert

# Extract required data from metamist
QUERY_PROJECT_ASSAYS = gql(
    """
    query ProjectAssays($datasetName: String!) {
    project(name: $datasetName) {
        samples {
            externalId
            meta
            assays {
                id
                meta
            }
        }
        }
    }
    """,
)


def query_metamist():
    """
    Query metamist for bioheart and tob-wgs datasets and map all tube ids to samples
    via call to create_default_dict.
    Later calls append_dictionaries(), which collates dicts from the two projects.

    :return: Dictionary of lists, mapping tube IDs to samples retrieved from metamist
            for both bioheart and tob-wgs datasets
    :rtype: defaultdict
    """
    # gives us json output from graphql query
    # Results are seen in graphql interface as a dictionary
    query_result_tob = query(QUERY_PROJECT_ASSAYS, variables={'datasetName': 'tob-wgs'})

    query_result_bioheart = query(
        QUERY_PROJECT_ASSAYS,
        variables={'datasetName': 'bioheart'},
    )

    # list of dictionaries of samples per External ID
    tob_samples = query_result_tob['project']['samples']
    bioheart_samples = query_result_bioheart['project']['samples']

    # Create default dicts for the two datasets: maps tube IDs to samples
    tob_kccg_id = 'KCCG FluidX tube ID'
    fluidx_to_assay_ids_tob = create_fluidx_to_assay_ids_dict(tob_samples, tob_kccg_id)
    bioheart_kccg_id = 'fluid_x_tube_id'
    fluidx_to_assay_ids_bioheart = create_fluidx_to_assay_ids_dict(
        bioheart_samples,
        bioheart_kccg_id,
    )

    return append_dictionaries(fluidx_to_assay_ids_tob, fluidx_to_assay_ids_bioheart)


def create_fluidx_to_assay_ids_dict(samples: list, kccg_id: str) -> defaultdict:
    """
    Creates defaultdict mapping kccg tube ids to all related assays.
    Flexible and can be applied to tob-wgs and bioheart.

    :return: defaultdict mapping kccg tubes to assay ids
    :rtype: defaultdict(list)
    """
    assays = []
    for sample in samples:
        assays.extend(sample['assays'])

    # Use default dict to group assays with fluid X tube metadata
    tube_to_assays_defaultdict = defaultdict(list)
    for assay in assays:
        # per assay, extract assay ID and fluid X tubeID
        # fluid tube id is the key; list contains assay IDs
        fluidx_tube_id = assay['meta'].get(kccg_id)
        if kccg_id.endswith('tube_id'):
            tube_to_assays_defaultdict[fluidx_tube_id.split('_')[1]].append(
                assay.get('id'),
            )
        else:
            tube_to_assays_defaultdict[fluidx_tube_id].append(assay['id'])
    return tube_to_assays_defaultdict


def append_dictionaries(
    fluidx_to_assay_ids_tob: defaultdict,
    fluidx_to_assay_ids_bioheart: defaultdict,
) -> defaultdict:
    """
    Append two defaultDict objects into a single dictionary.
    Verifies that there are no intersecting tube ids

    :return: When no intersection of keys
    """
    result = defaultdict(list)
    dicts_to_merge = [fluidx_to_assay_ids_tob, fluidx_to_assay_ids_bioheart]
    # First check to ensure that there are no overlapping keys
    set_tob = set(fluidx_to_assay_ids_tob.keys())
    set_bioheart = set(fluidx_to_assay_ids_bioheart.keys())

    if len(set_tob.intersection(set_bioheart)) == 0:
        # Append bioheart dict to tob dict
        for d in dicts_to_merge:
            for key, value in d.items():
                result[key].append(value)
    else:
        print(
            f'Error: these keys intersect: {set_tob.intersection(set_bioheart)}',
        )
        result = defaultdict()

    return result


def extract_excel():
    """
    Reads in metadata stored in excel spreadsheets. Appends these into
    a single dataframe and performs required transformations to columns
    to extract the desired data.

    :return: Fluid tube (key) IDS mapped to sequencing dates (value)
    :rtype: dict
    """
    workbook_names = [
        'scripts/metadata_enrichment/1K1K Sequencing Dates.xlsx',
        'scripts/metadata_enrichment/BioHEART Sequencing Dates.xlsx',
    ]
    df_list = []

    # Amalgamate data in all sheets listed in tob_sheet_names and bioheart_sheet_names
    for workbook in workbook_names:
        df_sheet_index = pd.read_excel(workbook, sheet_name=None, header=0)
        temp_df = pd.concat(df_sheet_index)
        df_list.append(temp_df)

    aggregated_df = pd.concat(df_list, ignore_index=True)

    # Separate accession date and fluidX_id in Sample Identifier column
    aggregated_df['fluidX_id'] = aggregated_df['Sample Identifier'].apply(
        lambda x: x.split('_')[1],
    )

    # Convert date to format YYYY-MM-DD
    aggregated_df['Date (YY/MM/DD)'] = pd.to_datetime(
        aggregated_df['Date (YY/MM/DD)'],
        format='%y%m%d',
    )
    aggregated_df['Date (YY/MM/DD)'] = aggregated_df['Date (YY/MM/DD)'].dt.strftime(
        '%Y-%m-%d',
    )

    # Insert sequencing date value into new dictionary fluidX_to_sequencing_date. Key is FluidX ID
    return pd.Series(
        aggregated_df['Date (YY/MM/DD)'].values,
        index=aggregated_df.fluidX_id,
    ).to_dict()


async def upsert_sequencing_dates(
    fluidx_to_assay_ids: defaultdict,
    fluidx_to_sequencing_date: dict,
):
    """
    Upserts meta {'fluid_x_tube_sequencing_date: x'} for each fluidx tube in metamist.
    Calls AssayUpsert endpoint update_assay_async with list of upserts: manages bulk calls asynchronously.

    :return: Results from bulk API assay upsert calls
    :rtype: Future
    """
    assay_api = AssayApi()
    assays_to_update = []

    for fluidx_id, assay_ids in fluidx_to_assay_ids.items():
        # Get sequencing date for each fluidX id from fluidX_to_sequencing_date dict
        if not fluidx_id:
            print(f'*****\nNO TUBE ID FOR : {fluidx_id} and assays {assay_ids}')
            continue

        try:
            # Create object required to update assay
            date = fluidx_to_sequencing_date[fluidx_id]
        except KeyError:
            print(f'Tube ID {fluidx_id} is not in the sequencing manifest')
        else:
            for assay_id in assay_ids:
                updated_assay = AssayUpsert(
                    id=assay_id,
                    meta={'fluid_x_tube_sequencing_date': date},
                )
                assays_to_update.append(
                    assay_api.update_assay_async(assay_upsert=updated_assay),
                )

    return await asyncio.gather(*assays_to_update)


def compare_tubes_metamist_excel(
    fluidx_to_assay_ids: defaultdict,
    fluidx_to_sequencing_date: dict,
):
    """
    Includes diagnostic code used to identify tubes missing between metamist and the
    provided sequencing date manifests. Prints output to console.
    """
    # Create a set from the fluidX identifiers only that were extracted from metamist
    metamist_set_fluidx = set(fluidx_to_assay_ids.keys())

    # Create a second set with fluidX identifiers extracted from excel files
    excel_set_fluidx = set(fluidx_to_sequencing_date.keys())

    # Find intersection/difference between two sets
    print(f'Count fluidX in Metamist {len(metamist_set_fluidx)}')
    print(f'Count fluidX in Excel {len(excel_set_fluidx)}')

    diff_metamist_excel = metamist_set_fluidx.difference(excel_set_fluidx)
    print(f'Diff metamist excel {diff_metamist_excel}')
    print(f'Len diff metamist excel {len(diff_metamist_excel)}')
    diff_excel_metamist = excel_set_fluidx.difference(metamist_set_fluidx)
    print(f'Diff excel metamist {diff_excel_metamist}')
    print(f'Len diff excel metamist {len(diff_excel_metamist)}')

    # Which assays are associated with the 'None' tube?
    assays_no_tubes = fluidx_to_assay_ids.get(None, 'empty')
    print(assays_no_tubes)
    print(len(assays_no_tubes[0]))


async def main():
    fluidx_to_assay_ids = query_metamist()
    if len(fluidx_to_assay_ids) != 0:
        fluidx_to_sequencing_date = extract_excel()
        # Uncomment this to compare the data in metamist to the manifests
        # compare_tubes_metamist_excel(fluidx_to_assay_ids, fluidx_to_sequencing_date) #noqa: ERA001
        await upsert_sequencing_dates(
            fluidx_to_assay_ids,
            fluidx_to_sequencing_date,
        )
    else:
        print('You have no FluidX_ids/assays to upsert')


if __name__ == '__main__':
    asyncio.new_event_loop().run_until_complete(main())
