#!/usr/bin/env python3

import hail as hl
from hail.methods import export_plink
from cpg_utils.hail_batch import dataset_path, init_batch, reference_path

# object containing variants within a 50K window on either side of the IGLL5 gene
MT = dataset_path('v0/IGLL5_50K_window.mt')

def main():
    # read and densify object
    init_batch()
    mt = hl.read_matrix_table(MT)
    mt = hl.experimental.densify(mt)

    # filter out low quality variants
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)

    # consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    # annotate using VEP
    vep_ht = reference_path('tob_wgs_vep/104/vep104.3_GRCh38.ht')
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    rv_mt = mt.filter_rows(mt.variant_qc.AF[1] < 0.05)

    # filter variants found to have regulatory effects
    filtered_mt = rv_mt.filtered_mt(hl.len(rv_mt.vep.regulatory_feature_consequences['biotype']) > 0)

    # export MT object to PLINK
    export_plink(filtered_mt, 'plink_files/igll5_rare_regulatory', ind_id=mt.s)

if __name__ == '__main__':
    main()