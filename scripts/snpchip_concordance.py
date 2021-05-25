"""Concordance between SNPchip genotypes and WGS variant calls"""

import os
import hail as hl

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-analysis/snpchip')

INPUT_VCF = 'gs://cpg-tob-wgs-analysis/snpchip/snpchip_vcf_rehead.vcf.bgz'
OUTPUT_MT = f'{output}/snpchip_hg38.mt'

hl.init()

mt = hl.import_vcf(INPUT_VCF)

# Explore MT
mt.count()
mt.rows().show()
mt.s.show(5)
mt.entry.show(5)


def liftover(x):
    """
    Liftover matrix table x from GRCh37 to GRCh38
    """
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover(
        'gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38
    )
    x = x.annotate_rows(new_locus=hl.liftover(x.locus, 'GRCh38'))
    x = x.filter_rows(hl.is_defined(x.new_locus))
    x = x.key_rows_by(locus=x.new_locus)
    return x


# Liftover to GRCh38
mt = liftover(mt)
mt.rows().show()

# Write lifted over MT to analysis bucket
mt.write(OUTPUT_MT)
