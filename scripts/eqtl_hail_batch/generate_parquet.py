"""Read in hail tables and export to parquet"""

import hail as hl
import click


@click.command()
@click.option('--input-path', help='A path prefix of where to input files are located')
def query(input_path):
    """Export tables to parquet"""

    hl.init(default_reference='GRCh38')

    # read in table and export as 
    table = input_path('eqtl_effect.ht')
    # write to parquet
    filename = input_path('eqtl_effect.parquet')
    table.to_spark().write.mode('append').parquet(filename)
