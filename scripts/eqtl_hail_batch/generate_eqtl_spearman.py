#!/usr/bin/env python3

"""Run Spearman rank correlation on SNPs and expression residuals"""

import os

import hail as hl
import hailtop.batch as hb
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from patsy import dmatrices  # pylint: disable=no-name-in-module
from scipy.stats import spearmanr
from cpg_utils.hail import copy_common_env, init_batch, remote_tmpdir
import click

DEFAULT_DRIVER_MEMORY = '4G'
DEFAULT_DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d2a9c316d6d752edb27623542c8a062db4466842-hail-0.2.73.devc6f6f09cec08'  # noqa: E501; pylint: disable=line-too-long
DRIVER_IMAGE = os.getenv('DRIVER_IMAGE', DEFAULT_DRIVER_IMAGE)

# TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v7.mt/'


def get_number_of_scatters(expression_df, geneloc_df):
    """get index of total number of genes"""

    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = list(expression_df.columns.values)[1:]
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(gene_ids)]

    return len(geneloc_df.index)


def filter_lowly_expressed_genes(expression_df):
    """Remove genes with low expression in all samples"""

    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
    # Number of individuals with non-zero expression
    genes_not_equal_zero = expression_df.iloc[:, 1:].values != 0
    n_expr_over_zero = pd.DataFrame(genes_not_equal_zero.sum(axis=0))
    percent_expr_over_zero = (n_expr_over_zero / len(expression_df.index)) * 100
    percent_expr_over_zero.index = expression_df.columns[1:]

    # Filter genes with less than 10 percent individuals with non-zero expression
    atleast10percent = percent_expr_over_zero[(percent_expr_over_zero > 10)[0]]
    sample_ids = expression_df['sampleid']
    expression_df = expression_df[atleast10percent.index]
    expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    return expression_df


def get_log_expression(expression_df):

    """get logged expression values"""

    expression_df = filter_lowly_expressed_genes(expression_df)

    # Prepare variables
    sample_ids = expression_df.iloc[:, 0]

    # log expression values
    to_log = expression_df.iloc[:, 1:].columns
    log_expression_df = expression_df[to_log].applymap(lambda x: np.log(x + 1))
    log_expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    return log_expression_df


def calculate_log_cpm(expression_df, output_prefix):
    """Calculate log cpm for each cell type and chromosome"""

    expression_df = filter_lowly_expressed_genes(expression_df)
    sample_ids = expression_df.iloc[:, 0]
    # remove sampleid column and get log expression
    # this can only be done on integers
    expression_df = expression_df.iloc[:, 1:]
    cpm_df = expression_df.apply(lambda x: (x / sum(x)) * 1000000, axis=0)
    log_cpm = np.log(cpm_df + 1)
    # add sampleids back in
    log_cpm = log_cpm.assign(sampleid=list(sample_ids))
    log_cpm.to_csv(os.path.join(output_prefix, f'log_cpm.tsv'), index=False)


def calculate_residuals(expression_df, covariate_df, output_prefix):
    """Calculate residuals for each gene in scatter"""

    log_expression_df = get_log_expression(expression_df)
    # Prepare variables
    gene_ids = list(log_expression_df.columns.values)[1:]
    sample_ids = log_expression_df.iloc[:, 0]

    # Calculate expression residuals
    def calculate_gene_residual(gene_id):
        """Calculate gene residuals"""
        gene = gene_id
        exprs_val = log_expression_df[['sampleid', gene]]
        test_df = exprs_val.merge(covariate_df, on='sampleid', how='left')
        test_df = test_df.rename(columns={test_df.columns[1]: 'expression'})
        y, x = dmatrices(
            'expression ~ sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + age + pf1 + pf2',
            test_df,
        )
        model = sm.OLS(y, x)
        residuals = list(model.fit().resid)
        return residuals

    residual_df = pd.DataFrame(list(map(calculate_gene_residual, gene_ids))).T
    residual_df.columns = gene_ids
    residual_df = residual_df.assign(sampleid=list(sample_ids))
    residual_df.to_csv(os.path.join(output_prefix, f'log_residuals.tsv'), index=False)

    return residual_df


# Run Spearman rank in parallel by sending genes in a batches
def run_spearman_correlation_scatter(
    idx,
    expression_df,
    genotype_df,
    geneloc_df,
    snploc_df,
    residuals_df,
    sampleid_keys,
):  # pylint: disable=too-many-locals
    """Run genes in scatter"""

    log_expression_df = get_log_expression(expression_df)

    # Prepare variables used to calculate Spearman's correlation
    gene_ids = list(log_expression_df.columns.values)[1:]
    # change genotype df sampleids from CPG internal IDs to OneK1K IDs
    onek1k_id = pd.merge(
        pd.DataFrame(genotype_df.sampleid),
        sampleid_keys,
        how='left',
        left_on='sampleid',
        right_on='InternalID',
    ).OneK1K_ID
    genotype_df.sampleid = onek1k_id
    # only keep samples that have rna-seq expression data
    genotype_df = genotype_df[genotype_df.sampleid.isin(log_expression_df.sampleid)]
    # FIXME: Only keep SNPs where 90% of individuals have values
    min_count = int(len(genotype_df.index) * 0.90)
    genotype_df = genotype_df.dropna(axis='columns', thresh=min_count)
    # filter snploc file to have the same snps as the genotype_df
    snploc_df = snploc_df[snploc_df.snpid.isin(genotype_df.columns)]

    # Get 1Mb sliding window around each gene
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(gene_ids)]
    geneloc_df = geneloc_df.assign(left=geneloc_df.start - 1000000)
    geneloc_df = geneloc_df.assign(right=geneloc_df.end + 1000000)

    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_symbol = df.gene_symbol
        gene_id = df.gene_id
        snp = df.snpid
        res_val = residuals_df[['sampleid', gene_symbol]]
        genotype_val = genotype_df[['sampleid', snp]]
        test_df = res_val.merge(genotype_val, on='sampleid', how='left')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        # set spearmanr calculation to perform the calculation ignoring nan values
        # this should be removed after resolving why NA values are in the genotype file
        coef, p = spearmanr(test_df['SNP'], test_df['residual'], nan_policy='omit')
        return (gene_symbol, gene_id, snp, coef, p)

    gene_info = geneloc_df.iloc[idx]
    print(gene_info)
    snps_within_region = snploc_df[
        snploc_df['pos'].between(gene_info['left'], gene_info['right'])
    ]
    gene_snp_df = snploc_df.merge(pd.DataFrame(snps_within_region))
    gene_snp_df = gene_snp_df.assign(
        gene_id=gene_info.gene_id, gene_symbol=gene_info.gene_name
    )
    spearman_df = pd.DataFrame(list(gene_snp_df.apply(spearman_correlation, axis=1)))
    spearman_df.columns = [
        'gene_symbol',
        'gene_id',
        'snpid',
        'spearmans_rho',
        'p_value',
    ]
    # add in global position and round
    locus = spearman_df.snpid.str.split(':', expand=True)[[0, 1]].agg(':'.join, axis=1)
    chrom = locus.str.split(':', expand=True)[0]
    bp = locus.str.split(':', expand=True)[1]
    spearman_df['locus'], spearman_df['chrom'], spearman_df['bp'] = [
        locus,
        chrom,
        bp,
    ]
    spearman_df['round'] = 1
    init_batch()
    t = hl.Table.from_pandas(spearman_df)
    t = t.annotate(global_bp=hl.locus(t.chrom, hl.int32(t.bp)).global_position())
    t = t.annotate(locus=hl.locus(t.chrom, hl.int32(t.bp)))
    # get alleles
    # mt = hl.read_matrix_table(TOB_WGS).key_rows_by('locus')
    t = t.key_by('locus')
    #    t = t.annotate(
    # alleles=mt.rows()[t.locus].alleles,
    # a1=mt.rows()[t.locus].alleles[0],
    # a2=hl.if_else(
    #     hl.len(mt.rows()[t.locus].alleles) == 2, mt.rows()[t.locus].alleles[1], 'NA'
    # ),
    #    )
    t = t.annotate(
        id=hl.str(':').join(
            [
                hl.str(t.chrom),
                hl.str(t.bp),
                # t.a1,
                # t.a2,
                t.gene_symbol,
                # result.db_key, # cell_type_id (eg nk, mononc)
                hl.str(t.round),
            ]
        )
    )
    spearman_df = t.to_pandas()
    return spearman_df


def merge_df_and_convert_to_string(*df_list):
    """Merge all Spearman dfs and convert to string using .to_string() on df"""
    merged_df: pd.DataFrame = pd.concat(df_list)
    pvalues = merged_df['p_value']
    fdr_values = pd.DataFrame(list(multi.fdrcorrection(pvalues))).iloc[1]
    merged_df = merged_df.assign(fdr=fdr_values)
    merged_df['fdr'] = merged_df.fdr.astype(float)
    return merged_df.to_string()


# Create click command line to enter dependency files
@click.command()
@click.option(
    '--expression', required=True, help='A sample x gene TSV of expression values'
)
@click.option('--genotype', required=True, help='A TSV of genotypes for each sample')
@click.option(
    '--geneloc', required=True, help='A TSV of start and end positions for each gene'
)
@click.option(
    '--snploc',
    required=True,
    help='A TSV of snp IDs with chromsome and position values for each',
)
@click.option(
    '--covariates', required=True, help='A TSV of covariates to calculate residuals'
)
@click.option(
    '--keys',
    required=True,
    help='A TSV of sample ids to convert external to internal IDs',
)  # pylint: disable=too-many-locals
@click.option(
    '--output-prefix',
    required=True,
    help='A path prefix of where to output files, eg: gs://MyBucket/output-folder/',
)
def main(
    expression: str,
    genotype,
    geneloc,
    snploc,
    covariates,
    keys,
    output_prefix: str,
):
    """
    Creates a Hail Batch pipeline for calculating EQTLs
    """
    dataset = os.getenv('CPG_DATASET')
    backend = hb.ServiceBackend(billing_project=dataset, remote_tmpdir=remote_tmpdir())
    batch = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

    # load in files literally to do the get_number of scatters
    expression_df_literal = pd.read_csv(expression, sep='\t')
    geneloc_df_literal = pd.read_csv(geneloc, sep='\t')

    # load files into a python job to avoid memory issues during a submission
    load_expression = batch.new_python_job('load-expression')
    load_expression.cpu(2)
    load_expression.memory('8Gi')
    load_expression.storage('2Gi')
    expression_df = load_expression.call(pd.read_csv, expression, sep='\t')
    load_genotype = batch.new_python_job('load-genotype')
    load_genotype.cpu(2)
    load_genotype.memory('8Gi')
    load_genotype.storage('2Gi')
    genotype_df = load_genotype.call(pd.read_csv, genotype, sep='\t')
    load_small_files = batch.new_python_job('load-small-files')
    load_small_files.memory('8Gi')
    load_small_files.storage('2Gi')
    geneloc_df = load_small_files.call(pd.read_csv, geneloc, sep='\t')
    snploc_df = load_small_files.call(pd.read_csv, snploc, sep='\t')
    covariate_df = load_small_files.call(pd.read_csv, covariates, sep=',')
    sampleid_keys = load_small_files.call(pd.read_csv, keys, sep='\t')
    calculate_residuals_job = batch.new_python_job('calculate-residuals')
    residuals_df = calculate_residuals_job.call(
        calculate_residuals,
        expression_df=expression_df,
        covariate_df=covariate_df,
        output_prefix=output_prefix,
    )
    calculate_log_cpm_job = batch.new_python_job('calculate-log-cpm')
    calculate_log_cpm_job.call(
        calculate_log_cpm,
        expression_df=expression_df,
        output_prefix=output_prefix,
    )

    spearman_dfs_from_scatter = []
    for idx in range(get_number_of_scatters(expression_df_literal, geneloc_df_literal)):
        # for i in range(5):
        j = batch.new_python_job(name=f'process_{idx}')
        j.cpu(2)
        j.memory('8Gi')
        j.storage('2Gi')
        copy_common_env(j)
        result: hb.resource.PythonResult = j.call(
            run_spearman_correlation_scatter,
            idx=idx,
            expression_df=expression_df,
            genotype_df=genotype_df,
            geneloc_df=geneloc_df,
            snploc_df=snploc_df,
            residuals_df=residuals_df,
            sampleid_keys=sampleid_keys,
        )
        spearman_dfs_from_scatter.append(result)

    merge_job = batch.new_python_job(name='merge_scatters')
    merge_job.cpu(2)
    merge_job.memory('8Gi')
    merge_job.storage('2Gi')
    result_second = merge_job.call(
        merge_df_and_convert_to_string, *spearman_dfs_from_scatter
    )
    corr_result_output_path = os.path.join(output_prefix + 'correlation_results.csv')
    batch.write_output(result_second.as_str(), corr_result_output_path)
    batch.run(wait=False)


if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    main()
