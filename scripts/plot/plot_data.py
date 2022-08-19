#!/usr/bin/env python3


"""
Test plotting
"""

import hail as hl
from cpg_utils.hail_batch import output_path
from bokeh.io.export import get_screenshot_as_png


def main():
    """
    Run vep using main.py wrapper
    """

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table('gs://cpg-tob-wgs-test/tob_wgs_vep/v1/vep105_GRCh38.mt')
    mt = hl.filter_intervals(
        mt,
        [hl.parse_locus_interval('chr22:23704425-23802743', reference_genome='GRCh38')],
    )
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    p1 = hl.plot.histogram(mt.variant_qc.AF[1])
    p1_filename = output_path(f'histogram_maf_post_filter.png', 'web')
    with hl.hadoop_open(p1_filename, 'wb') as f:
        get_screenshot_as_png(p1).save(f, format='PNG')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
