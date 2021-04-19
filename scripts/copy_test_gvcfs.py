"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')

SUBSET_RE = ''

for batch in [0, 1]:
    subprocess.run(
        [
            'gsutil',
            'cp',
            f'gs://cpg-tob-wgs-main/gvcf/batch{batch}/{SUBSET_RE}?.g.vcf.gz*',
            f'gs://cpg-tob-wgs-test/gvcf/batch{batch}/',
        ],
        check=False,
    )

# Copy metadata subset.
subprocess.run(
    f'gsutil cat gs://cpg-tob-wgs-main/batch{batch}/*.csv | '
    f'grep -e "sample\\|{SUBSET_RE}" | '  # Also copy header.
    f'gsutil cp - gs://cpg-tob-wgs-main/gvcf/batch{batch}/',
    shell=True,
    check=False,
)
