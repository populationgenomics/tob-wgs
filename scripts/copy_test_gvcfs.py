"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')

for bn in [0, 1]:
    subprocess.run(
        [
            'gsutil',
            'cp',
            f'gs://cpg-tob-wgs-main/gvcf/batch{bn}/*.g.vcf.gz*',
            f'gs://cpg-tob-wgs-test/gvcf/batch{bn}/',
        ],
        check=False,
    )

    # Copy metadata subset.
    subprocess.run(
        [
            'gsutil',
            'cp',
            f'gs://cpg-tob-wgs-main/gvcf/batch{bn}/*.csv',
            f'gs://cpg-tob-wgs-test/gvcf/batch{bn}/',
        ],
        check=False,
    )
