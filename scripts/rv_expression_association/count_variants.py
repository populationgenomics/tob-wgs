#!/usr/bin/env python3  

# This script aims to count the total number of variants (non-ref)
# from the whole TOB-WGS dataset that are biallelic SNVs
# and rare (MAF < 5%)

import hail as hl
from cpg_utils.hail_batch import dataset_path, init_batch

# full TOB-WGS object
MT = dataset_path('mt/v7.mt')


def main():  
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(MT)
    # densify
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)                                   # remove hom-ref
        & (mt.n_unsplit_alleles == 2)                                 # biallelic
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))                   # SNVs
    )

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0))

    print(mt.count())

    
if __name__ == '__main__':
    main()
    
