#!/usr/bin/env python3  

import click
import hail as hl
from cpg_utils.hail_batch import dataset_path, init_batch


@click.command()
@click.option('--mt-path', required=True)  # 'mt/v7.mt'
@click.option('--genes', required=True)    # 'LMNA'
def subset_variants(
    mt_path: str,
    genes: str,
):  

    # read ht objects to extract variant names
    for gene in genes:
        ht_object_filename = f'gs://cpg-tob-wgs-main-analysis/tob_wgs_rv/pseudobulk_rv_association/summary_hts/{gene}_rare_promoter_summary.ht'
        ht = hl.read_table(ht_object_filename)
        mt

    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(dataset_path(mt_path))
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
    subset_variants()