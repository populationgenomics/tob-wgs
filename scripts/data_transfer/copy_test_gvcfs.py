"""Copies a subset of gVCFs from the main bucket to the test bucket."""

import os
import subprocess

QC_METADATA_FILE_NAME = 'R_210315_BINKAN1_1K1KDNA_M001.csv'

SAMPLES = [
    'TOB1520',
    'TOB1521',
    'TOB1522',
    'TOB1523',
    'TOB1524',
    'TOB1525',
    'TOB1526',
    'TOB1527',
    'TOB1528',
    'TOB1529',
    'TOB1530',
    'TOB1531',
    'TOB1532',
    'TOB1533',
    'TOB1534',
    'TOB1535',
    'TOB1536',
    'TOB1537',
    'TOB1538',
    'TOB1640',
]

# for sn in SAMPLES:
#     subprocess.run(
#         [
#             'gsutil',
#             'cp',
#             f'gs://cpg-tob-wgs-main/gvcf/batch1/{sn}.g.vcf.gz*',
#             'gs://cpg-tob-wgs-test/gvcf/batch1/',
#         ],
#         check=False,
#     )

# Copy metadata subset.
subprocess.run(
    f'gsutil cat gs://cpg-tob-wgs-main-metadata/batch1/{QC_METADATA_FILE_NAME} | '
    f'grep -P "sample|{"|".join(SAMPLES)}" | '  # Also copy header.
    f'gsutil cp - gs://cpg-tob-wgs-test-metadata/batch1/{QC_METADATA_FILE_NAME}',
    shell=True,
    check=False,
)
