#!/usr/bin/env python3  
# pylint: disable=missing-module-docstring,missing-function-docstring,no-value-for-parameter,import-error,no-name-in-module

# This script aims to count the total number of variants (non-ref)
# from the whole TOB-WGS dataset that are biallelic SNVs
# and rare (MAF < 5%)

import click
import hail as hl
from cpg_utils.hail_batch import dataset_path, init_batch


@click.command()
@click.option('--mt-path', required=True)  # 'mt/v7.mt'
def count_variants(
    mt_path: str,
):  
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
    count_variants()
    
