#!/usr/bin/env python3


"""
Run VEP on the hail mt
"""


import click
import hail as hl
from cpg_utils.hail_batch import output_path


@click.command()
@click.option('--mt', required=True, help='Hail matrix table to run VEP on')
def main(mt: str):
    """
    Run vep using main.py wrapper
    """

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(mt)
    ht = mt.rows()
    # filter to biallelic loci only
    ht = ht.filter(hl.len(ht.alleles) == 2)
    ht = ht.filter(ht.alleles[1] != '*')
    vep = hl.vep(mt, config='file:///vep_data/vep-gcloud.json')
    vep_path = output_path('vep104_GRCh38.mt')
    vep.write(vep_path)


if __name__ == '__main__':
    main()
