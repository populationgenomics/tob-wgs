from metamist.graphql import query, gql
from collections import defaultdict
import pandas as pd
from metamist.apis import AssayApi
from metamist.models import AssayUpsert
import asyncio
import csv

# Extract required data from metamist
QUERY_TOB = gql(
    """
    query enrich_tob_seq_date {
    project(name: "tob-wgs") {
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
    """)

QUERY_BIOHEART = gql(
    """
    query enrich_bioheart_seq_date {
    project(name: "bioheart") {
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
    """)

QUERY_SEQ_GR_ACTIVE = gql(
    """
    query check_cram_assays {
        project(name: "tob-wgs") {
            participants {
            externalId
            samples {
                sequencingGroups(activeOnly: {contains: true}) {
                id
                assays {
                    id
                    meta
                }
                }
            }
            }
        }
    }
    """
    )

def query_metamist():
    # gives us json output from graphql query
    # Results are seen in graphql interface as a dictionary
    query_result = query(QUERY_TOB)
    query_result_bioheart = query(QUERY_BIOHEART)

    # list of dictionaries of samples per External ID
    tob_samples = query_result['project']['samples']
    bioheart_samples = query_result_bioheart['project']['samples']

    # TOB-WGS
    # Create list of all assays from tob_samples (at an individual level)
    assays = []
    for sample in tob_samples:
        assays.extend(sample.get('assays'))

    # Use default dict to group assays with fluid X tube metadata
    fluidX_to_assay_ids = defaultdict(list)
    for assay in assays:
        # per assay, extract assay ID and fluid X tubeID
        # fluid tube id is the key; list contains assay IDs
        fluidX_tube_id = assay.get('meta').get('KCCG FluidX tube ID')
        # Check that FluidX tube is not null TODO
        fluidX_to_assay_ids[fluidX_tube_id].append(assay.get('id'))

    # BIOHEART
    assays_bioheart = []
    for sample in bioheart_samples:
        assays_bioheart.extend(sample.get('assays'))

    # Use default dict to group assays with fluid X tube metadata
    fluidX_to_assay_ids_bioheart = defaultdict(list)
    for assay in assays_bioheart:
        # per assay, extract assay ID and fluid X tubeID
        # fluid tube id is the key; list contains assay IDs
        fluidX_tube_id_bioheart = assay.get('meta').get('fluid_x_tube_id')
        fluidX_to_assay_ids_bioheart[fluidX_tube_id_bioheart].append(assay.get('id'))

    return fluidX_to_assay_ids, fluidX_to_assay_ids_bioheart

def extract_excel():
    tob_workbook_names = ['scripts/metadata_enrichment/1K1K Sequencing Dates.xlsx', 'scripts/metadata_enrichment/BioHEART Sequencing Dates.xlsx']
    # aggregated_df = pd.DataFrame()
    df_list = list()

    # Amalgamate data in all sheets listed in tob_sheet_names and bioheart_sheet_names
    for workbook in tob_workbook_names:
        df_sheet_index = pd.read_excel(workbook, sheet_name=None, header=0)
        temp_df = pd.concat(df_sheet_index)
        df_list.append(temp_df)

    aggregated_df = pd.concat(df_list, ignore_index=True)
    # Separate accession date and fluidX_id in Sample Identifier column
    aggregated_df['accession_date'] = aggregated_df['Sample Identifier'].apply(lambda x: x.split('_')[0])
    aggregated_df['fluidX_id'] = aggregated_df['Sample Identifier'].apply(lambda x: x.split('_')[1])

    # Insert accession value into new dictionary fluidX_to_sequencing_date. Key is FluidX ID
    fluidX_to_sequencing_date = pd.Series(aggregated_df.accession_date.values,index=aggregated_df.fluidX_id).to_dict()
    return fluidX_to_sequencing_date

def upsert_sequencing_dates(fluidX_to_assay_ids, fluidX_to_sequencing_date):
    # construct API update calls
    # Iterate through the fluidX_to_assay_IDs dictionary because this is representative of what's already in metamist
    # That is: fluidX_to_assay_IDs groups assays by fluidX_tube_id

    assay_API = AssayApi()
    api_calls_to_gather = []

    for fluidX_id, assay_ids in fluidX_to_assay_ids.items(): 
    # Get sequencing date for each fluidX id from fluidX_to_sequencing_date dict
        if fluidX_id: 
            try:
                date = fluidX_to_sequencing_date[fluidX_id]
            except KeyError:
                print(f'Tube ID {fluidX_id} is not in provided sequencing manifest')
            else:
                for assay_id in assay_ids:
                    # Create API call and save to 
                    api_calls_to_gather.append(assay_API.update_assay_async(AssayUpsert(id=assay_id, meta={'sequencing_date': date})))
                    # assay_API.update_assay(AssayUpsert(id=assay_id, meta={'sequencing_date': date}))

        else:
            print(f'*****\nNO TUBE ID FOR : {fluidX_id} and assays {assay_ids}')

    # TODO: confirm asyncio documentation; create async function 
    # return await asyncio.gather(api_calls_to_gather)

def compare_tubes_metamist_excel(fluidX_to_assay_ids, fluidX_to_sequencing_date): 
    # Create a set from the fluidX identifiers only that were extracted from metamist
    metamist_set_fluidX = set(fluidX_to_assay_ids.keys())

    # Create a second set with fluidX identifiers extracted from excel files
    excel_set_fluidX = set(fluidX_to_sequencing_date.keys())

    # Find intersection/difference between two sets
    print(f'Count fluidX in Metamist {len(metamist_set_fluidX)}')
    print(f'Count fluidX in Excel {len(excel_set_fluidX)}')

    print(f'Diff metamist excel {metamist_set_fluidX.difference(excel_set_fluidX)}')
    diff_excel_metamist = excel_set_fluidX.difference(metamist_set_fluidX)
    print(f'Diff excel metamist {diff_excel_metamist}')
    print(f'Len diff excel metamist {len(diff_excel_metamist)}')
    # Save diff_excel_metamist to csv

    with open('./scripts/metadata_enrichment/bioheart_tubes_missing_in_Metamist.csv', 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(['tube_id'])
        for item in diff_excel_metamist:
            writer.writerow([item])


def is_active_sequencing_group():
    query_result = query(QUERY_SEQ_GR_ACTIVE)

    active_samples = query_result['project']['participants']
    external_ids_to_active_assay_ids = defaultdict(list)

    participants = []
    # group all active samples by external ID
    for participant in active_samples:
        assay_ids = []
        # participants.extend(participant.get('samples'))
        external_id = participant.get('externalId')
        # Get the id for all our active sequencing groups
    #     assay_ids = participant.get('samples').get('sequencingGroups').get('assays').get('id')
        assay_ids = participant.get('samples').get('sequencingGroups')
        print(assay_ids)
    #     for id in assay_ids:
    #         # First, you need a list of all the 
    #             external_ids_to_active_assay_ids[external_id].append(id)

                
    # for i in external_ids_to_active_assay_ids:
    #     print(i)

    # Create list of all assays from tob_samples (at an individual level)
    # assays = []
    # for sample in tob_samples:
    #     assays.extend(sample.get('assays'))

    # # Use default dict to group assays with fluid X tube metadata
    # fluidX_to_assay_ids = defaultdict(list)
    # for assay in assays:
    #     # per assay, extract assay ID and fluid X tubeID
    #     # fluid tube id is the key; list contains assay IDs
    #     fluidX_tube_id = assay.get('meta').get('KCCG FluidX tube ID')
    #     fluidX_to_assay_ids[fluidX_tube_id].append(assay.get('id'))



if __name__ == '__main__':
    fluidX_to_assay_ids, fluidX_tube_id_bioheart = query_metamist()
    fluidX_to_sequencing_date = extract_excel()
    upsert_sequencing_dates(fluidX_to_assay_ids, fluidX_to_sequencing_date)

    # Exploration only 
    # compare_tubes_metamist_excel(fluidX_to_assay_ids, fluidX_to_sequencing_date)
    # is_active_sequencing_group()