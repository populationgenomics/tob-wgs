# We need metamist query and parsing for excel file

from metamist.graphql import query, gql
from collections import defaultdict
import pandas as pd
from metamist.apis import AssayApi
from metamist.models import AssayUpsert

# Get what we want out of metamist
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

# gives us json output as seen in graphql interface as dictionary
query_result = query(QUERY)

# list full of dictionaries of samples 
tob_samples = query_result['project']['samples']

assays = []
for sample in tob_samples:
   assays.extend(sample.get('assays'))

fluidX_to_assay_IDs = defaultdict(list)
for assay in assays:
   # extract assay ID and fluid X tubeID
   # fluid tube id is the key; list contains assay IDs
   fluidX_tube_id = assay.get('meta').get('KCCG FluidX tube ID')
   fluidX_to_assay_IDs[fluidX_tube_id].append(assay.get('id'))

tob_sheet_names = [] # Make list

# For all sheets listed in tob_sheet_names and bioheart_sheet_names
for sheet in tob_sheet_names:
   df_sheet_index = pd.read_excel('sample.xlsx', sheet_name=sheet)

   # Start with blank dataframe and stack each sheet onto the end: make a giant dataframe
   # for column 1, rip off the accession date (first 7 chars)
   # Insert this value into a new dictionary fluidX_to_sequencing_date. Values are date values (whichever is easiest)

# construct API update calls
# Iterate through the fluidX_to_assay_IDs dictionary because this is what's already in metamist

assay_API = AssayApi()

for fluidX_id, assay_IDs in fluidX_to_assay_IDs.items(): 
   # Get sequencing date for each fluidX id from our dictionary that doesn't exist yet
    date = 00000000
    for assay_id in assay_IDs:
        # comment this out until you're ready to execute the script
        # assay_API.update_assay(AssayUpsert(id=assay_id, meta={'sequencing_date': date}))
        print("Validate content with print statements or move this into a notebook")

# https://pythonbasics.org/read-excel/
# https://github.com/populationgenomics/metamist/blob/dev/tob_metadata_harmonisation.py
# update_assayAsync can let us bundle everything together. Develop script in 'snail way' and then adapt for batch/async implementation. 
