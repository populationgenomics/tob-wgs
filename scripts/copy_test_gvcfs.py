"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')

subprocess.run(
    ['gsutil', 'cp', 'gs://cpg-tob-wgs-main/v0/TOB15[2-3]?.g.vcf.gz*', output]
)
