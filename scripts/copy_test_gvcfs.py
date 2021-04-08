"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')

SUBSET_RE = 'TOB15[2-3]'
QC_METADATA_FILE_NAME = 'R_210315_BINKAN1_1K1KDNA_M001.csv'

subprocess.run(
    ['gsutil', 'cp', f'gs://cpg-tob-wgs-main/v0/{SUBSET_RE}?.g.vcf.gz*', output]
)

# Copy metadata subset.
subprocess.run(
    f'gsutil cat gs://cpg-tob-wgs-main/v0/{QC_METADATA_FILE_NAME} | '
    f'grep -e "sample\\|{SUBSET_RE}" | '  # Also copy header.
    f'gsutil cp - {output}/{QC_METADATA_FILE_NAME}',
    shell=True,
)
