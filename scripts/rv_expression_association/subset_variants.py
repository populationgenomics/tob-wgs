#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,missing-function-docstring,missing-module-docstring

import click

import hail as hl

from cpg_utils.hail_batch import dataset_path, init_batch, output_path


@click.command()
@click.option('--input-mt-path', required=True)  # 'mt/v7.mt'
@click.option(
    '--genes',
    required=True,
    help='List of genes to consider. Space separated, as one argument.',
)  # 'LMNA'
@click.option('--output-mt-prefix', required=True)  # 'significant_gene_burden_max/'
def subset_variants(
    input_mt_path: str,
    genes: str,
    output_mt_prefix: str,
):
    init_batch()

    # read hail matrix table object (WGS data)
    mt = hl.read_matrix_table(dataset_path(input_mt_path))
    mt = hl.experimental.densify(mt)

    # read ht objects and subset mt to specific variants
    genes_of_interest = genes.split(' ')
    for gene in genes_of_interest:
        ht_object_filename = dataset_path(
            f'tob_wgs_rv/pseudobulk_rv_association/summary_hts/{gene}_rare_promoter_summary.ht',
            'analysis',
        )
        ht = hl.read_table(ht_object_filename)
        # merge objects
        mt_tmp = mt.semi_join_rows(ht)
        # save mt
        output_mt_name = f'{output_mt_prefix}{gene}.mt'
        output_mt_path = output_path(output_mt_name, 'analysis')
        mt_tmp.write(output_mt_path, overwrite=True)


if __name__ == '__main__':
    subset_variants()
