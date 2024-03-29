#!/usr/bin/env python3
# ruff: noqa: PLR2004
"""Plot call rate distribution for TOB-WGS dataset"""

from bokeh.io.export import get_screenshot_as_png

import hail as hl

from cpg_utils.hail_batch import dataset_path, init_batch, output_path

MT = dataset_path('mt/v7.mt/')


def main():
    """plot call rate distribution"""

    init_batch()
    mt = hl.read_matrix_table(MT)
    samples = mt.s.collect()
    n_samples = len(samples)
    mt = hl.experimental.densify(mt)
    # filter to biallelic loci only
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    # filter out variants that didn't pass the VQSR filter
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    # Calculate call rate and add to mt
    mt = mt.annotate_rows(call_rate=hl.agg.sum(hl.is_defined(mt.GT)) / n_samples)
    # plot distribution
    p1 = hl.plot.histogram(
        mt.call_rate,
        range=(0, 1),
        bins=30,
        legend='Call Rate',
        title='Call Rate Distribution',
    )
    p1_filename = output_path('histogram_call_rate.png', 'web')
    with hl.hadoop_open(p1_filename, 'wb') as f:
        get_screenshot_as_png(p1).save(f, format='PNG')
    p2 = hl.plot.cumulative_histogram(
        mt.call_rate,
        range=(0, 1),
        bins=30,
        legend='Call Rate',
        title='Call Rate Distribution',
    )
    p2_filename = output_path('cumulative_histogram_call_rate.png', 'web')
    with hl.hadoop_open(p2_filename, 'wb') as f:
        get_screenshot_as_png(p2).save(f, format='PNG')


if __name__ == '__main__':
    main()
