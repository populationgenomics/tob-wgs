# flake8: noqa: PD015
"""
Create metadata sheet for TOB samples
Run using the command:
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "scrna-seq/grch38_association_files/metadata" \
--description "TOB metadata" python3 create_metadata_sheet.py
"""

import logging

import pandas as pd

from cpg_utils.hail_batch import output_path

METADATA = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/metadata/metadata.csv'
)
OUTLIERS = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/metadata/outliers.csv'
)
SAMPLEID_KEYS = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
)

metadata = pd.read_csv(METADATA)
outliers = pd.read_csv(OUTLIERS)
combined_metadata = metadata.merge(
    pd.DataFrame(outliers),
    how='left',
    left_on='CPG_ID',
    right_on='samples',
)
combined_metadata = combined_metadata.rename(columns={'samples': 'population'})
# make sure all 51 outliers properly mapped to a CPG ID
n_outliers = len(combined_metadata.population.dropna())
logging.info(f'{n_outliers} samples are population outliers')
combined_metadata.population = combined_metadata.population.fillna('NFE')
combined_metadata.loc[
    combined_metadata['population'].str.contains('CPG'),
    'population',
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
merged_cpg_onek1k.to_csv(output_path('keys_metadata_sheet.csv'), index=False)
