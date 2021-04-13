"""Moves gVCFs from the upload bucket to the main bucket.
Eventually this will be automated by the upload processor."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-main/gvcf/batch')

subprocess.run(
    [
        'gsutil',
        'mv',
        'gs://cpg-tob-wgs-upload/*.g.vcf.gz*',
        'gs://cpg-tob-wgs-upload/*.csv',  # QC metadata.
        output,
    ]
)
