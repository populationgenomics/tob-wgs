#!/usr/bin/env python3

"""Generate genotype dfs for the association analysis"""

import hail as hl

TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v7.mt/'

hl.init(default_reference='GRCh38')

mt = hl.read_matrix_table(TOB_WGS)
# filter out variants that didn't pass the VQSR filter
mt = mt.filter_rows(hl.is_missing(mt.filters))
# VQSR does not filter out low quality genotypes. Filter these out
mt = mt.filter_entries(mt.GQ <= 20, keep = False)
# filter out samples with a genotype call rate > 0.8 
n_samples = mt.count_cols()
call_rate = 0.8
mt = mt.filter_rows(hl.agg.sum(hl.is_missing(mt.GT)) > (n_samples * call_rate), keep=False)
# filter out variants with MAF < 0.01
mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
# turn files into output that we need
# turn into plink, split up by chromosome


