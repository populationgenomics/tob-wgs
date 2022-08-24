#!/usr/bin/env python3

"""Entry point for the analysis runner."""

from bokeh.io.export import get_screenshot_as_png

# from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, init_batch, output_path, reference_path
import hail as hl

VEP_MT = dataset_path('tob_wgs_vep/v1/vep105_GRCh38.mt')


def main():
    init_batch()

    mt = hl.read_matrix_table(VEP_MT)
    # mt = hl.experimental.densify(mt)
    mt = hl.filter_intervals(
        mt,
        [hl.parse_locus_interval('chr22:23704425-23802743', reference_genome='GRCh38')],
    )
    mt = hl.experimental.densify(mt)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    p1 = hl.plot.histogram(mt.variant_qc.AF[1])
    p1_filename = output_path('histogram_maf_pre_filter.png', 'web')
    with hl.hadoop_open(p1_filename, 'wb') as f:
        get_screenshot_as_png(p1).save(f, format='PNG')


    sample = 'CPG18'
    donor_mt = mt.filter_cols(mt.s == sample)
    donor_mt = donor_mt.filter_rows(hl.agg.any(donor_mt.GT.is_non_ref()))
    mt = mt.semi_join_rows(donor_mt.rows())
    p2 = hl.plot.histogram(mt.variant_qc.AF[1])
    p2_filename = output_path('histogram_maf_post_filter_old.png', 'web')
    with hl.hadoop_open(p2_filename, 'wb') as f:
        get_screenshot_as_png(p2).save(f, format='PNG')

    ref_ht = reference_path('seqr/v0-1/combined_reference_data_grch38-2.0.4.ht')
    donor_mt = donor_mt.annotate_rows(
        gnomad_genomes=ref_ht[donor_mt.row_key].gnomad_genomes,
    )

    p3 = hl.plot.histogram(donor_mt.gnomad_genomes.AF_POPMAX_OR_GLOBAL)
    p3_filename = output_path('histogram_maf_post_filter_gnomad.png', 'web')
    with hl.hadoop_open(p3_filename, 'wb') as f:
        get_screenshot_as_png(p3).save(f, format='PNG')


if __name__ == '__main__':
    main()
