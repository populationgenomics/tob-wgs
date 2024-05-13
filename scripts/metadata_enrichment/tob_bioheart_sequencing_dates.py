import csv
from collections import defaultdict

import pandas as pd

from metamist.apis import AssayApi
from metamist.graphql import gql, query

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


# TODO: Create separate function for default dict creation
def query_metamist():
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

    # TOB-WGS
    # Create list of all assays from tob_samples (at an individual level)
    assays = []
    for sample in tob_samples:
        assays.extend(sample['assays'])

    # Use default dict to group assays with fluid X tube metadata
    fluidx_to_assay_ids_tob = defaultdict(list)
    for assay in assays:
        # per assay, extract assay ID and fluid X tubeID
        # fluid tube id is the key; list contains assay IDs
        fluidx_tube_id = assay['meta'].get('KCCG FluidX tube ID')
        # Check that FluidX tube is not null TODO
        fluidx_to_assay_ids_tob[fluidx_tube_id].append(assay['id'])

    # BIOHEART
    assays_bioheart = []
    for sample in bioheart_samples:
        assays_bioheart.extend(sample.get('assays'))

    # Use default dict to group assays with fluid X tube metadata
    fluidx_to_assay_ids_bioheart = defaultdict(list)
    for assay in assays_bioheart:
        # per assay, extract assay ID and fluid X tubeID
        # fluid tube id is the key; list contains assay IDs
        fluidx_tube_id_bioheart = assay.get('meta').get('fluid_x_tube_id')
        fluidx_to_assay_ids_bioheart[fluidx_tube_id_bioheart.split('_')[1]].append(
            assay.get('id'),
        )

    return append_dictionaries(fluidx_to_assay_ids_tob, fluidx_to_assay_ids_bioheart)


def append_dictionaries(
    fluidx_to_assay_ids_tob: defaultdict,
    fluidx_to_assay_ids_bioheart: defaultdict,
) -> defaultdict:
    """
    Append two defaultDict objects into a single dictionary.
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
        return result
    else:  # noqa: RET505
        print(
            f'Error: these keys intersect: {set_tob.intersection(set_bioheart)}',
        )
        return defaultdict()


# TODO: combine pd.apply() into a single expression
def extract_excel():
    tob_workbook_names = [
        'scripts/metadata_enrichment/1K1K Sequencing Dates.xlsx',
        'scripts/metadata_enrichment/BioHEART Sequencing Dates.xlsx',
    ]
    df_list = []

    # Amalgamate data in all sheets listed in tob_sheet_names and bioheart_sheet_names
    for workbook in tob_workbook_names:
        df_sheet_index = pd.read_excel(workbook, sheet_name=None, header=0)
        temp_df = pd.concat(df_sheet_index)
        df_list.append(temp_df)

    aggregated_df = pd.concat(df_list, ignore_index=True)
    # Separate accession date and fluidX_id in Sample Identifier column

    # TODO: Can you combine these into a single expression?
    aggregated_df['accession_date'] = aggregated_df['Sample Identifier'].apply(
        lambda x: x.split('_')[0],
    )
    aggregated_df['fluidX_id'] = aggregated_df['Sample Identifier'].apply(
        lambda x: x.split('_')[1],
    )

    # Insert accession value into new dictionary fluidX_to_sequencing_date. Key is FluidX ID
    return pd.Series(
        aggregated_df.accession_date.values,
        index=aggregated_df.fluidX_id,
    ).to_dict()


# TODO: combine api calls and
def upsert_sequencing_dates(
    fluidx_to_assay_ids: defaultdict,
    fluidx_to_sequencing_date: dict,
):
    # construct API update calls
    # Iterate through the fluidX_to_assay_IDs dictionary because this is representative of what's already in metamist
    # That is: fluidX_to_assay_IDs groups assays by fluidX_tube_id

    assay_api = AssayApi()  # noqa: F841
    # api_calls_to_gather = [] #noqa: ERA001

    for fluidx_id, assay_ids in fluidx_to_assay_ids.items():
        # Get sequencing date for each fluidX id from fluidX_to_sequencing_date dict
        if fluidx_id:
            try:
                date = fluidx_to_sequencing_date[fluidx_id]  # noqa: F841
            except KeyError:
                print(f'Tube ID {fluidx_id} is not in the sequencing manifest')
            else:
                for assay_id in assay_ids:
                    print(f'API call added for {assay_id}')
                    # api_calls_to_gather.append(assay_API.update_assay_async(AssayUpsert(id=assay_id, meta={'sequencing_date': date}))) # noqa: ERA001
                    # assay_API.update_assay(AssayUpsert(id=assay_id, meta={'sequencing_date': date})) # noqa: ERA001

        else:
            print(f'*****\nNO TUBE ID FOR : {fluidx_id} and assays {assay_ids}')
    # TODO: confirm asyncio documentation; create async function
    # return await asyncio.gather(api_calls_to_gather) #noqa: ERA001


def compare_tubes_metamist_excel(
    fluidx_to_assay_ids: defaultdict,
    fluidx_to_sequencing_date: dict,
):
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

    # Sort diff_metamist_excel (determine if chronological)
    print(sorted(filter(None, diff_metamist_excel)))

    # Save diff_metamist_excel to csv
    with open(
        './scripts/metadata_enrichment/tubes_missing_in_manifests.csv',
        'w',
    ) as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(['tube_id'])
        for item in diff_metamist_excel:
            writer.writerow([item])


if __name__ == '__main__':
    fluidx_to_assay_ids = query_metamist()
    if len(fluidx_to_assay_ids) != 0:
        fluidx_to_sequencing_date = extract_excel()
        # Exploration only
        # compare_tubes_metamist_excel(fluidx_to_assay_ids, fluidx_to_sequencing_date) # noqa: ERA001
    else:
        print('You have no FluidX_ids/assays to upsert')
