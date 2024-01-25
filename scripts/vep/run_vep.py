#!/usr/bin/env python3
# flake8: noqa: PLR2004

"""
Run VEP on the hail mt
"""


import click

import hail as hl

from cpg_utils.hail_batch import dataset_path, output_path


@click.command()
@click.option('--mt', 'mt_path', required=True, help='Hail matrix table to run VEP on')
@click.option('--vep-version', help='Version of VEP', default='104.3')
def main(mt_path: str, vep_version: str):
    """
    Run vep using generate_genotype_df_batch.py wrapper
    """

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(dataset_path(mt_path))
    ht = mt.rows()
    # filter to biallelic loci only
    ht = ht.filter(hl.len(ht.alleles) == 2)
    # remove starred alleles, as this results in an error in VEP
    # see https://discuss.hail.is/t/vep-output-variant-not-found-in-original-variants/1148
    ht = ht.filter(ht.alleles[1] != '*')
    vep = hl.vep(ht, config='file:///vep_data/vep-gcloud.json')
    vep_path = output_path(f'vep{vep_version}_GRCh38.ht')
    vep.write(vep_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
