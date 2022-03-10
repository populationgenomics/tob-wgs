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
    locus_alleles = mt.rows().select().to_pandas(flatten=True)
    print(locus_alleles.head())
    # expand locus to two columns and rename 
    locus_alleles = locus_alleles.astype(str)
    print(locus_alleles.head())
    x = locus_alleles['locus'].str.split(':', expand=True)
    print(x.head())
    print(pd.concat([x, locus_alleles['alleles']], axis=1).head())
    locus_alleles = pd.concat([locus_alleles['locus'].str.split(':', expand=True), locus_alleles['alleles']], axis=1)
    print(locus_alleles.head())
    locus_alleles.columns = ['locus.contig', 'locus.position', 'alleles']
    print(locus_alleles.head())
    # save each chromosome to an individual file
    for chr in set(locus_alleles['locus.contig']): 
        print(chr)
        locus_alleles.loc[locus_alleles['locus.contig'] == chr].to_parquet(f'gs://cpg-tob-wgs-test/kat/v0/{chr}_test.parquet')

if __name__ == '__main__':
    query()
