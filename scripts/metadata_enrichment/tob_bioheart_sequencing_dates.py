from metamist.graphql import query, gql
from collections import defaultdict
import pandas as pd
from metamist.apis import AssayApi
from metamist.models import AssayUpsert

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
fluidX_to_assay_IDs = defaultdict(list)
for assay in assays:
   # per assay, extract assay ID and fluid X tubeID
   # fluid tube id is the key; list contains assay IDs
   fluidX_tube_id = assay.get('meta').get('KCCG FluidX tube ID')
   fluidX_to_assay_IDs[fluidX_tube_id].append(assay.get('id'))

tob_workbook_names = ['scripts/metadata_enrichment/1K1K Sequencing Dates.xlsx', 'scripts/metadata_enrichment/BioHEART Sequencing Dates.xlsx']
# aggregated_df = pd.DataFrame()
df_list = list()

# Amalgamate data in all sheets listed in tob_sheet_names and bioheart_sheet_names
for workbook in tob_workbook_names:
    df_sheet_index = pd.read_excel(workbook, sheet_name=None, header=0)
    temp_df = pd.concat(df_sheet_index)
    df_list.append(temp_df)

aggregated_df = pd.concat(df_list, ignore_index=True)
# Separate accession date and fluidX_id
aggregated_df['accession_date'] = aggregated_df['Sample Identifier'].apply(lambda x: x.split('_')[0])
aggregated_df['fluidX_id'] = aggregated_df['Sample Identifier'].apply(lambda x: x.split('_')[1])

# Insert accession value into new dictionary fluidX_to_sequencing_date. Key is FluidX ID
fluidX_to_sequencing_date = pd.Series(aggregated_df.accession_date.values,index=aggregated_df.fluidX_id).to_dict()

# construct API update calls
# Iterate through the fluidX_to_assay_IDs dictionary because this is representative of what's already in metamist
# That is: fluidX_to_assay_IDs groups assays by fluidX_tube_id

assay_API = AssayApi()

for fluidX_id, assay_ids in fluidX_to_assay_IDs.items(): 
    print(f'F ID: {fluidX_id}')
   # Get sequencing date for each fluidX id from fluidX_to_sequencing_date dict
    date = fluidX_to_sequencing_date[fluidX_id]
    print(date)
    for assay_id in assay_ids:
        # comment this out until you're ready to execute the script
        # assay_API.update_assay(AssayUpsert(id=assay_id, meta={'sequencing_date': date}))
        print(f'F ID: {fluidX_id} assay ID: {assay_id}')

# https://pythonbasics.org/read-excel/
# https://github.com/populationgenomics/metamist/blob/dev/tob_metadata_harmonisation.py
# update_assayAsync can let us bundle everything together. Develop script in 'snail way' and then adapt for batch/async implementation. 
