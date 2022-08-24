#!/usr/bin/env python3

"""Entry point for the analysis runner."""

from bokeh.io.export import get_screenshot_as_png

# from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, init_batch, output_path
import hail as hl

VEP_MT = dataset_path('tob_wgs_vep/v1/vep105_GRCh38.mt')


def main():
    init_batch()

    mt = hl.read_matrix_table(VEP_MT)
    mt = head(mt, 10000)
    mt = hl.experimental.densify(mt)
    mt = hl.variant_qc(mt)
    mt = hl.filter_intervals(
        mt,
        [hl.parse_locus_interval('chr22:23704425-23802743', reference_genome='GRCh38')],
    )
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    p1 = hl.plot.histogram(mt.variant_qc.AF[1])
    p1_filename = output_path('histogram_maf_pre_filter.png', 'web')
    with hl.hadoop_open(p1_filename, 'wb') as f:
        get_screenshot_as_png(p1).save(f, format='PNG')

    # p2_filename = output_path('histogram_plot2.png', 'web')
    # with hl.hadoop_open(p2_filename, 'wb') as f:
    #     get_screenshot_as_png(p1).save(f, format='PNG')

    sample = 'CPG18'
    donor_mt = mt.filter_cols(mt.s == sample)
    donor_mt = donor_mt.filter_rows(hl.agg.any(donor_mt.GT.is_non_ref()))
    mt = mt.semi_join_rows(donor_mt.rows())
    p2 = hl.plot.histogram(mt.variant_qc.AF[1])
    p2_filename = output_path('histogram_maf_post_filter.png', 'web')
    with hl.hadoop_open(p2_filename, 'wb') as f:
        get_screenshot_as_png(p2).save(f, format='PNG')


if __name__ == '__main__':
    main()
