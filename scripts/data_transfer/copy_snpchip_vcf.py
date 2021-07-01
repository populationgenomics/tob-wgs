#!/usr/bin/env python3

"""
Copies raw SNPchip genotype data tarball and Plink-generated VCF
(transferred from GWCCG in 2021-Mar-17) to the tob-wgs main bucket.
"""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('snpchip/v1')

INPUT_BUCKET = 'gs://cpg-tob-wgs-snpchipdata'
INPUT_VCF = f'{INPUT_BUCKET}/data/onek1k_pre_imputation_genotypes.vcf.gz'
INPUT_TARBALL = f'{INPUT_BUCKET}/data/OneK1K_scRNA_GenotypingData.tar.gz'

subprocess.run(
    ['gsutil', 'cp', INPUT_VCF, INPUT_TARBALL, "gs://cpg-tob-wgs-main/snpchip/v1/"],
    check=False,
)
