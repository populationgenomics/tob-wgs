#!/usr/bin/env python3

"""
Copies GRCh38 SNPchip MatrixTable to the tob-wgs main bucket.
"""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('snpchip/v1')

INPUT_MT = 'gs://cpg-tob-wgs-test/snpchip/v1/snpchip_grch38.mt'

subprocess.run(['gsutil', '-m', 'cp', '-r', INPUT_MT, output], check=False)
