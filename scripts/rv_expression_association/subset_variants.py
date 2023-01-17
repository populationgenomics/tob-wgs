#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,missing-function-docstring

import click
import hail as hl
from cpg_utils.hail_batch import dataset_path, init_batch


@click.command()
@click.option('--input-mt-path', required=True)  # 'mt/v7.mt'
@click.option('--genes', required=True)  # 'LMNA'
@click.option(
    '--output-mt-path', required=True
)  # 'gs://cpg-tob-wgs-main-analysis/tob_wgs_rv/pseudobulk_rv_association/significant_genes_burden_max.mt'
def subset_variants(
    input_mt_path: str,
    genes: str,
    output_mt_path: str,
):

    # read ht objects to extract variant names
    variants = list()
    genes_of_interest = genes.split(' ')
    for gene in genes_of_interest:
        ht_object_filename = f'gs://cpg-tob-wgs-main-analysis/tob_wgs_rv/pseudobulk_rv_association/summary_hts/{gene}_rare_promoter_summary.ht'
        ht = hl.read_table(ht_object_filename)
        variants.append(ht.locus.collect())  # I am nearly sure the syntax is wrong here

    # read hail matrix table object (WGS data)
    init_batch()  # does this need to be above?
    mt = hl.read_matrix_table(dataset_path(input_mt_path))
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(mt.locus in variants)

    # save mt
    mt.write(
        output_mt_path, overwrite=True
    )  # how do I specify the analysis bucket again?


if __name__ == '__main__':
    subset_variants()

