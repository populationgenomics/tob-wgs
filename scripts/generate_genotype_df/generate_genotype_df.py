#!/usr/bin/env python3

"""Generate genotype dfs for the association analysis"""

import hail as hl
import os
import hailtop.batch as hb
from cpg_utils.hail import dataset_path, output_path, copy_common_env, init_query_service

TOB_WGS = dataset_path('mt/v7.mt/')


def generate_genotypes():
    """Generate genotype dfs for each chromosome"""
    init_query_service()

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
    # select alleles and locus (contig and position must be selected separately),
    # then convert to pandas
    t = mt.rows()
    t = t.key_by(contig=t.locus.contig, position=t.locus.position)
    t = t.select(t.alleles)
    n_rows = t.count()
    print(f'n rows = {n_rows}')
    pd = t.to_pandas(flatten=True)
    # expand locus to two columns and rename 
    # save each chromosome to an individual file
    for chr in set(pd['contig']): 
        pd.loc[pd['contig'] == chr].to_parquet(output_path(f'tob_genotype_maf01_{chr}.parquet'))


if __name__ == '__main__':
    dataset = os.getenv('CPG_DATASET')
    access_level = os.getenv('CPG_ACCESS_LEVEL')
    backend = hb.ServiceBackend(billing_project=dataset, remote_tmpdir=f'cpg-{dataset}-{access_level}')
    b = hb.Batch()
    j = b.new_python_job('generate-genotypes')
    j.image(os.getenv('CPG_DRIVER_IMAGE'))
    j.memory('32Gi')
    copy_common_env(j)
    j.call(generate_genotypes)

    b.run(wait=False)
