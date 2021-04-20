"""Moves gVCFs from the upload bucket to the main bucket.
Eventually this will be automated by the upload processor."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-main/gvcf/batch')

# Make sure this isn't adding to a previous batch.
assert (
    subprocess.run(['gsutil', 'ls', f'{output}/*.g.vcf.gz*'], check=False).returncode
    == 1
)

subprocess.run(
    [
        'gsutil',
        'mv',
        'gs://cpg-tob-wgs-upload/*.g.vcf.gz*',
        'gs://cpg-tob-wgs-upload/*.csv',  # QC metadata.
        output,
    ],
    check=True,
)
