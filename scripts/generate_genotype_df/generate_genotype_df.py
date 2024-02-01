#!/usr/bin/env python3
# ruff: noqa: PLR2004
"""Generate genotype dfs for the association analysis"""

import hail as hl

from cpg_utils.hail_batch import dataset_path, output_path

TOB_WGS = dataset_path('mt/v7.mt/')
FREQ_TABLE = dataset_path('joint-calling/v7/variant_qc/frequencies.ht/', 'analysis')


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
    mt = mt.filter_rows(
        hl.agg.sum(hl.is_missing(mt.GT)) > (n_samples * call_rate),
        keep=False,
    )
    # filter out variants with MAF < 0.01
    ht = hl.read_table(FREQ_TABLE)
    mt = mt.annotate_rows(freq=ht[mt.row_key].freq)
    mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
    # select alleles and locus (contig and position must be selected separately),
    # then convert to pandas
    t = mt.rows()
    t = t.key_by(contig=t.locus.contig, position=t.locus.position)
    t = t.select(t.alleles)
    pd = t.to_pandas(flatten=True)
    # expand locus to two columns and rename
    # save each chromosome to an individual file
    for chromosome in set(pd['contig']):
        pd.loc[pd['contig'] == chromosome].to_parquet(
            output_path(f'tob_genotype_maf01_{chromosome}.parquet'),
        )


if __name__ == '__main__':
    query()
