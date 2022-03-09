#!/usr/bin/env python3

"""Generate genotype dfs for the association analysis"""

import hail as hl
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v7.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(TOB_WGS)
    # filter out variants that didn't pass the VQSR filter
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)

    # VQSR does not filter out low quality genotypes. Filter these out
    mt = mt.filter_entries(mt.GQ <= 20, keep=False)
    # filter out samples with a genotype call rate > 0.8 (as in the gnomAD supplementary paper)
    n_samples = mt.count_cols()
    call_rate = 0.8
    mt = mt.filter_rows(hl.agg.sum(hl.is_missing(mt.GT)) > (n_samples * call_rate), keep=False)
    # filter out variants with MAF < 0.01
    mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
    tob_wgs_path = output_path('tob_wgs_maf01.parquet')
    mt.rows().to_pandas()[['locus.contig', 'locus.position', 'alleles']].to_parquet(tob_wgs_path)


if __name__ == '__main__':
    query()
