#!/usr/bin/env python3

"""Entry point for the analysis runner."""

from bokeh.io.export import get_screenshot_as_png

# from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, init_batch, output_path
# , reference_path
import hail as hl

# VEP_MT = dataset_path('tob_wgs_vep/v1/vep105_GRCh38.mt')
VEP_MT = dataset_path('mt/v7.mt')


def main():
    init_batch()

    mt = hl.read_matrix_table(VEP_MT)
    mt = hl.filter_intervals(
        mt,
        [hl.parse_locus_interval('chr22:23000000-24000000', reference_genome='GRCh38')],
    )
    mt = hl.experimental.densify(mt)

    mt = hl.variant_qc(mt)
    p1 = hl.plot.histogram(mt.variant_qc.call_rate, legend='Call Rate')
    p1_filename = output_path('histogram_variant_call_rate.png', 'web')
    with hl.hadoop_open(p1_filename, 'wb') as f:
        get_screenshot_as_png(p1).save(f, format='PNG')

    mt = hl.sample_qc(mt)
    p2 = hl.plot.histogram(mt.sample_qc.call_rate, legend='Call Rate')
    p2_filename = output_path('histogram_sample_call_rate.png', 'web')
    with hl.hadoop_open(p2_filename, 'wb') as f:
        get_screenshot_as_png(p2).save(f, format='PNG')
    
    p3 = hl.plot.histogram(mt.sample_qc.gq_stats.mean, legend='Mean Sample GQ')
    p3_filename = output_path('histogram_sample_gq.png', 'web')
    with hl.hadoop_open(p3_filename, 'wb') as f:
        get_screenshot_as_png(p3).save(f, format='PNG')


if __name__ == '__main__':
    main()
