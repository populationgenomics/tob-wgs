#!/usr/bin/env python3

import hail as hl
from bokeh.io import show
from bokeh.io.export import get_screenshot_as_png
from analysis_runner import output_path

mt = hl.read_matrix_table('gs://cpg-tob-wgs-test/tob_wgs_vep/v1/vep105_GRCh38.mt')
mt = hl.experimental.densify(mt)
mt = hl.variant_qc(mt)
mt = hl.filter_intervals(
    mt, [hl.parse_locus_interval('chr22:23704425-23802743', reference_genome='GRCh38')]
)
mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
p1 = hl.plot.histogram(mt.variant_qc.AF[1])
# show(p1)
p1_filename = output_path(f'histogram_maf_pre_filer.png', 'web')
with hl.hadoop_open(p1_filename, 'wb') as f:
    get_screenshot_as_png(p1).save(f, format='PNG')

sample = 'CPG18'
donor_mt = mt.filter_cols(mt.s == sample)
donor_mt = donor_mt.filter_rows(hl.agg.any(donor_mt.GT.is_non_ref()))
mt = mt.semi_join_rows(donor_mt.rows())
p2 = hl.plot.histogram(mt.variant_qc.AF[1])
# show(p2)
p2_filename = output_path(f'histogram_maf_post_filer.png', 'web')
with hl.hadoop_open(p2_filename, 'wb') as f:
    get_screenshot_as_png(p2).save(f, format='PNG')
