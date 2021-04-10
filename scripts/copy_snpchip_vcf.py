"""
Copies SNPchip VCF and scRNA tarball (transferred from GWCCG in 2021-Mar-17)
to the tob-wgs temporary bucket.
"""

import os
import subprocess

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-temporary/')

# INPUT_VCF = 'gs://cpg-tob-wgs-snpchipdata/data/onek1k_pre_imputation_genotypes.vcf.gz'
INPUT_TARBALL = 'gs://cpg-tob-wgs-snpchipdata/data/OneK1K_scRNA_GenotypingData.tar.gz'

subprocess.run(['gsutil', 'cp', INPUT_TARBALL, output], check=False)
