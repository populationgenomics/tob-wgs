"""
Script enriches existing assays in Metamist with sequencing dates
Sequencing dates are extracted from an external .csv file, which maps
dates to fluidx kccg tube ids.
"""

import asyncio
from collections import defaultdict

import click
import pandas as pd
from google.cloud import storage

from metamist.apis import AssayApi
from metamist.graphql import gql, query
from metamist.models import AssayUpsert
from metamist.parser.generic_metadata_parser import run_as_sync

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


def query_metamist(project: str):
    """
    Query metamist for project and map all tube ids to samples
    via call to create_default_dict.

    :param project str: Target project name. Originally provided as CLI arg
    :return: Dictionary of lists, mapping tube IDs to samples retrieved from metamist
    :rtype: defaultdict
    """
    # gives us json output from graphql query
    # Results are seen in graphql interface as a dictionary
    query_result = query(QUERY_PROJECT_ASSAYS, variables={'datasetName': project})

    # list of dictionaries of samples per External ID
    samples = query_result['project']['samples']

    # Create default dicts for the project: maps tube IDs to samples
    kccg_id = 'fluid_x_tube_id' if project == 'bioheart' else 'KCCG FluidX tube ID'

    return create_fluidx_to_assay_ids_dict(samples, kccg_id)


def create_fluidx_to_assay_ids_dict(samples: list, kccg_id: str) -> defaultdict:
    """
    Creates defaultdict mapping kccg tube ids to all related assays.
    Flexible and can be applied to tob-wgs and bioheart.

    :param samples list: List of all samples associated with defined project in Metamist
    :param kccg_id str: Tube IDs are not uniformly defined. kccg_id == id used in target project
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
        # Keys for bioheart/tob projects are formatted differently
        if kccg_id.endswith('tube_id'):
            tube_to_assays_defaultdict[fluidx_tube_id.split('_')[1]].append(
                assay.get('id'),
            )
        else:
            tube_to_assays_defaultdict[fluidx_tube_id].append(assay['id'])

    return tube_to_assays_defaultdict


def extract_excel(workbook_names: tuple):
    """
    Reads in metadata stored in excel spreadsheets. Appends these into
    a single dataframe and performs required transformations to columns
    to extract the desired data.

    :param workbook_names tuple: Name of files containing enrichment metadata,
        provided as CLI args. Initially there were multiple workbooks, hence
        use of 'for' loop.
    :return: Fluid tube (key) IDS mapped to sequencing dates (value)
    :rtype: dict
    """
    df_list = []
    storage_client = storage.Client()

    # Amalgamate data in all sheets listed in tob_sheet_names and bioheart_sheet_names
    for workbook in workbook_names:
        # Get workbook file information in GCP
        split_workbook_name = workbook.split('/')
        bucket_name = split_workbook_name[2]
        file_name = split_workbook_name[3]
        bucket = storage_client.get_bucket(bucket_name)
        blob = bucket.get_blob(file_name)
        workbook_file = blob.download_as_string()

        # Collate all workbooks as df into a single list for concatenation outside loop
        df_sheet_index = pd.read_excel(workbook_file, sheet_name=None, header=0)
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

    :param fluidx_to_assay_ids defaultdict: Tube IDs : sample ids in Metamist mapping. Stored as dictionary of
        lists with Tube ID as key.
    :param fluidx_to_sequencing_date dict: Tube:sequencing date mapping, as extracted from .xlsx metadata file.
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

    :param fluidx_to_assay_ids defaultdict: Tube IDs : sample ids in Metamist mapping. Stored as dictionary of
        lists with Tube ID as key.
    :param fluidx_to_sequencing_date dict: Tube:sequencing date mapping, as extracted from .xlsx metadata file.
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


@click.command()
@click.argument('project', nargs=1)
@click.argument('manifests', nargs=-1)
@click.option('--debug', '-d', is_flag=True)
@run_as_sync
async def main(project: str, manifests: tuple, debug: bool):
    """
    Sample CLI command:
        python scripts/metadata_enrichment/tob_bioheart_sequencing_dates.py <project/dataset name> \
        <PATH TO METADATA FILE IN GCP>

    :param project str: Target project name. Either bioheart or tob-wgs
    :param workbook_names tuple: Name of files containing enrichment metadata.
    :param debug bool: Boolean flag. True will run script in debug mode, which
        prints out summary data for project.
    """
    fluidx_to_assay_ids = query_metamist(project)
    if len(fluidx_to_assay_ids) != 0:
        fluidx_to_sequencing_date = extract_excel(manifests)
        if not debug:
            await upsert_sequencing_dates(
                fluidx_to_assay_ids,
                fluidx_to_sequencing_date,
            )
        else:
            compare_tubes_metamist_excel(
                fluidx_to_assay_ids,
                fluidx_to_sequencing_date,
            )
    else:
        print('You have no FluidX_ids/assays to upsert')


if __name__ == '__main__':
    asyncio.new_event_loop().run_until_complete(main())
