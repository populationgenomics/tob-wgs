"""Moves CRAMs from the upload bucket to the archive bucket.
Eventually this will be automated by the upload processor."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-archive/cram/batch')

# Make sure this isn't adding to a previous batch.
assert (
    subprocess.run(['gsutil', 'ls', f'{output}/*.cram*'], check=False).returncode == 1
)

subprocess.run(['gsutil', 'mv', 'gs://cpg-tob-wgs-upload/*.cram*', output], check=True)
