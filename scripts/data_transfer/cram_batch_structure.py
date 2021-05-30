"""Change the CRAM directory structure from versions to batches."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output == 'gs://cpg-tob-wgs-archive/cram/batch0/'

subprocess.run(['gsutil', 'mv', 'gs://cpg-tob-wgs-archive/v0/*', output], check=True)

# Copy missing sample TOB1521.
subprocess.run(
    ['gsutil', 'mv', 'gs://cpg-tob-wgs-upload/TOB1521.cram*', output], check=True
)
