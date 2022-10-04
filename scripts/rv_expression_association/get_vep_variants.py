#!/usr/bin/env python3  # pylint: disable=missing-module-docstring

import logging
import hail as hl
from hail.methods import export_plink
from cpg_utils.hail_batch import dataset_path, init_batch, output_path

# object containing variants within a 50K window on either side of the VPREB3 gene
MT = dataset_path('v0/VPREB3_50K_window.mt')
VEP_HT = dataset_path('v0/VPREB3_50K_window_vep.ht')

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


def main():  # pylint: disable=missing-function-docstring
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(MT)
    logging.info(f'Number of variants in window: {mt.count()[0]}')

    # densify
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)                                   # biallelic
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))                   # SNVs
    )

    # annotate using VEP
    vep_ht = hl.read_table(VEP_HT)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)
    mt_path = output_path('densified_qced_snps_only_vep_annotated.mt', 'tmp')
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f'Number of QC-passing, biallelic SNPs: {mt.count()[0]}')

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    rv_mt = mt.filter_rows((mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0))
    rv_mt_path = output_path('rare_variants.mt', 'tmp')
    rv_mt = rv_mt.checkpoint(rv_mt_path, overwrite=True)  # checkpoint
    logging.info(f'Number of rare variants (freq<5%): {rv_mt.count()[0]}')

    # filter variants found to have regulatory effects
    filtered_mt = rv_mt.filter_rows(
        hl.len(rv_mt.vep.regulatory_feature_consequences['biotype']) > 0
    )
    filtered_mt_path = output_path('regulatory_variants.mt', 'tmp')
    filtered_mt = filtered_mt.checkpoint(filtered_mt_path, overwrite=True)  # checkpoint
    logging.info(
        f'Number of rare variants (freq<5%) with regulatory consequences: {filtered_mt.count()[0]}'
    )

    # export this as a Hail table for downstream analysis
    ht_filename = output_path('vpreb3_rare_regulatory_summary.ht')
    ht = filtered_mt.rows()
    ht.write(ht_filename)  # will this work? should i select important colunmns only?

    # export MT object to PLINK (all regulatory variants)
    filtered_mt_path = output_path('plink_files/vpreb3_rare_regulatory')
    export_plink(filtered_mt, filtered_mt_path, ind_id=filtered_mt.s)

    # filter to promoter variants only
    filtered_mt2 = filtered_mt.filter_rows(
        filtered_mt.vep.regulatory_feature_consequences['biotype'][0] == 'promoter'
    )
    # export to PLINK
    filtered_mt_path2 = output_path('plink_files/vpreb3_rare_promoter')
    export_plink(filtered_mt2, filtered_mt_path2, ind_id=filtered_mt2.s)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
