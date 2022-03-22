"""Create metadata sheet for TOB samples"""

import pandas as pd

METADATA = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/metadata/metadata.csv'
)
OUTLIERS = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/metadata/outliers.csv'
)
SAMPLEID_KEYS = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
)
OUTPUT_DIR = 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/metadata'

metadata = pd.read_csv(METADATA)
outliers = pd.read_csv(OUTLIERS)
combined_metadata = pd.merge(
    metadata, pd.DataFrame(outliers), how='left', left_on='CPG_ID', right_on='samples'
)
combined_metadata = combined_metadata.rename(columns={'samples': 'population'})
# make sure all 51 outliers properly mapped to a CPG ID
len(combined_metadata.population.dropna())
# 51
combined_metadata.population = combined_metadata.population.fillna('NFE')
combined_metadata.loc[
    combined_metadata['population'].str.contains('CPG'), 'population'
] = 'outlier'
# add in oneK1K IDs
sampleid_keys = pd.read_csv(SAMPLEID_KEYS, sep='\t')
# Merge on TOB IDs, as you can get the OneK1K ID repeated (i.e., a CPG_ID that has two different
# idenentifiers, but the same TOB identifier, will have a repeated OneK1K ID)
merged_cpg_onek1k = pd.merge(
    sampleid_keys,
    combined_metadata,
    how='right',
    right_on='TOB_ID',
    left_on='ExternalID',
)
# save
output_path = f'{OUTPUT_DIR}/keys_metadata_sheet.csv'
merged_cpg_onek1k.to_csv(output_path, index=False)
