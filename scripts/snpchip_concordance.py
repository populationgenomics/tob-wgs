"""Concordance between SNPchip genotypes and WGS variant calls"""

import os
import subprocess
import hail as hl

# Rename sample columns that correspond to TOB IDs (51 columns do not have a TOB ID)
DATA_DIR = os.path.join('nogit/data/snpchip')
SNPCHIP_VCF_RAW = os.path.join(
    f'{DATA_DIR}/raw/onek1k_pre_imputation_genotypes_2021-Mar-17.vcf.bgz'
)
SNPCHIP_VCF_REHEAD = os.path.join(f'{DATA_DIR}/processed/hail/0-rehead.vcf.bgz')
SNPCHIP_ID2TOBID = os.path.join(
    f'{DATA_DIR}/processed/sample2tobid_bcftools_reheader_2021-May-17.tsv'
)

subprocess.run(
    [
        'bcftools',
        'reheader',
        '-s',
        SNPCHIP_ID2TOBID,
        SNPCHIP_VCF_RAW,
        '-o',
        SNPCHIP_VCF_REHEAD,
    ],
    check=False,
)

# Explore VCF in Hail
hl.init()

SNPCHIP_MT = os.path.join(f'{DATA_DIR}/processed/hail/0-rehead.mt')

hl.import_vcf(SNPCHIP_VCF_REHEAD).write(SNPCHIP_MT)
mt = hl.read_matrix_table(SNPCHIP_MT)

# Explore MT
mt.count()  # 483,482 x 1,034
mt.rows().show()  # show CHR:POS, [REF,ALT], ID, QUAL, FILTER, INFO
mt.s.show(5)  # show first 5 samples
mt.entry.show(5)  # show first 5 GTs


def liftover(x):
    """
    Liftover from GRCh37 to GRCh38
    """
    ref37 = hl.get_reference('GRCh37')
    ref38 = hl.get_reference('GRCh38')
    ref37.add_liftover(
        os.path.join(f'{DATA_DIR}/reference/grch37_to_grch38.over.chain.gz'), ref38
    )
    x = x.annotate_rows(new_locus=hl.liftover(x.locus, 'GRCh38'))
    x = x.filter_rows(hl.is_defined(x.new_locus))
    x = x.key_rows_by(locus=x.new_locus)
    return x


mt = liftover(mt)
mt.rows().show()
