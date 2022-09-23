#!/usr/bin/env python3

import hail as hl
import numpy as np
import pandas as pd
from hail.methods import export_plink
from cpg_utils.hail_batch import dataset_path, init_batch#, reference_path

# object containing variants within a 50K window on either side of the IGLL5 gene
MT = dataset_path('v0/IGLL5_50K_window.mt')
VEP_HT = dataset_path('v0/IGLL5_50K_window_vep.ht')

def main():
    # read and densify object
    init_batch()
    mt = hl.read_matrix_table(MT)
    mt = hl.experimental.densify(mt)

    print(mt.count())

    # filter out low quality variants
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)

    # consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    # annotate using VEP
    vep_ht = hl.read_table(VEP_HT)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)

    print(mt.count())

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    rv_mt = mt.filter_rows(mt.variant_qc.AF[1] < 0.05)
    rv_mt = rv_mt.filter_rows(rv_mt.variant_qc.AF[1] > 0)

    print(rv_mt.count())

    # filter variants found to have regulatory effects
    filtered_mt = rv_mt.filter_rows(hl.len(rv_mt.vep.regulatory_feature_consequences['biotype']) > 0)
    print(filtered_mt.count())

    filtered_rrv_mt = filtered_mt.filter_rows(filtered_mt.variant_qc.AF[1] < 0.01)
    print(filtered_rrv_mt.count())

    filtered_0maf_mt = filtered_mt.filter_rows(filtered_mt.variant_qc.AF[1] == 0)
    print(filtered_0maf_mt.count())

    # print out stats
    # types of regulatory variants
    biotypes = pd.Series(filtered_mt.vep.regulatory_feature_consequences['biotype'].collect())
    print(biotypes.value_counts())

    # MAF distribution
    mafs = pd.Series(filtered_mt.variant_qc.AF[1].collect())
    print(np.nanmin(mafs))
    print(np.nanmax(mafs))
    print(np.nanmean(mafs))

    # export MT object to PLINK
    export_plink(filtered_mt, 'plink_files/igll5_rare_regulatory', ind_id=mt.s)

    # generate SetID file

if __name__ == '__main__':
    main()