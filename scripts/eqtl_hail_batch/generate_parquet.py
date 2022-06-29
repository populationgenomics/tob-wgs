"""Read in hail tables and export to parquet"""

import hail as hl
from cpg_utils.hail_batch import output_path


def query():
    """Export tables to parquet"""

    hl.init(default_reference='GRCh38')

    # read in table and export as 
    table = output_path('gene_expression.ht')
    # write to parquet
    filename = output_path('eqtl_effect.parquet')
    table.to_spark().write.mode('append').parquet(filename)
