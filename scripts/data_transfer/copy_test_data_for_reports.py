"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')


subprocess.run(
    [
        'gsutil',
        'cp',
        'gs://cpg-tob-wgs-main/gvcf/batch4/R_210315_BINKAN1_1K1KDNA_M004.csv',
        'gs://cpg-tob-wgs-test/reporting/R_210315_BINKAN1_1K1KDNA_M004.csv',
    ],
    check=False,
)

subprocess.run(
    [
        'gsutil',
        'cp',
        'gs://cpg-tob-wgs-main/joint-calling/v2/meta.tsv',
        'gs://cpg-tob-wgs-test/reporting/meta.tsv',
    ],
    check=False,
)
