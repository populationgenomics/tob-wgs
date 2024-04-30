from metamist.graphql import query, gql
from collections import defaultdict
import pandas as pd
from metamist.apis import AssayApi
from metamist.models import AssayUpsert
import asyncio

# Extract required data from metamist
QUERY = gql(
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

def query_metamist():
    # gives us json output from graphql query 
    # Results are seen in graphql interface as a dictionary
    query_result = query(QUERY)

    # list of dictionaries of samples per External ID
    tob_samples = query_result['project']['samples']

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
        fluidX_to_assay_ids[fluidX_tube_id].append(assay.get('id'))
    
    return fluidX_to_assay_ids

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
    # await asyncio.gather(api_calls_to_gather)

if __name__ == '__main__':
    fluidX_to_assay_ids = query_metamist()
    fluidX_to_sequencing_date = extract_excel()
    upsert_sequencing_dates(fluidX_to_assay_ids, fluidX_to_sequencing_date)
