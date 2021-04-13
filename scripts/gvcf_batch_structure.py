"""Change the gVCF directory structure from versions to batches."""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output == 'gs://cpg-tob-wgs-main/gvcf/batch0/'

subprocess.run(['gsutil', 'mv', 'gs://cpg-tob-wgs-main/v0/*', output], check=True)

# Copy missing sample TOB521.
subprocess.run(
    ['gsutil', 'mv', 'gs://cpg-tob-wgs-upload/TOB1521.g.vcf.gz*', output], check=True
)
