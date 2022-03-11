#!/usr/bin/env python3

"""Generate genotype dfs for the association analysis"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v7.mt/')

def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(TOB_WGS)
    mt = hl.experimental.densify(mt)
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
    # select only locus and alleles, which are the keys, then convert to pandas
    mt = mt.annotate_rows(contig = mt.locus.contig, position=mt.locus.position).key_rows_by('alleles')
    locus_alleles = mt.rows().select('contig','position').to_pandas(flatten=True)
    # expand locus to two columns and rename 
    # save each chromosome to an individual file
    for chr in set(locus_alleles['contig']): 
        locus_alleles.loc[locus_alleles['contig'] == chr].to_parquet(output_path(f'tob_genotype_maf01_{chr}.parquet'))

if __name__ == '__main__':
    query()
