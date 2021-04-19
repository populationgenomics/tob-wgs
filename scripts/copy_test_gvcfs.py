"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')

SUBSET_RE = ''

subprocess.run(
    [
        'gsutil',
        'cp',
        f'gs://cpg-tob-wgs-main/gvcf/batch*/{SUBSET_RE}?.g.vcf.gz*',
        output,
    ],
    check=False,
)

# Copy metadata subset.
subprocess.run(
    f'gsutil cat gs://cpg-tob-wgs-main/batch1/*.csv | '
    f'grep -e "sample\\|{SUBSET_RE}" | '  # Also copy header.
    f'gsutil cp - {output}/',
    shell=True,
    check=False,
)
