#!/usr/bin/env python3

"""Run Spearman rank correlation on SNPs and expression residuals"""

import os

import hail as hl
import hailtop.batch as hb
import pandas as pd
import numpy as np
import statsmodels.api as sm
from patsy import dmatrices  # pylint: disable=no-name-in-module
from scipy.stats import spearmanr
from cpg_utils.hail_batch import (
    dataset_path,
    output_path,
    copy_common_env,
    init_batch,
    remote_tmpdir,
)
from cloudpathlib import AnyPath
import click
from multipy.fdr import qvalue

DEFAULT_DRIVER_MEMORY = '4G'
DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16'
assert DRIVER_IMAGE

TOB_WGS = dataset_path('mt/v7.mt/')
FREQ_TABLE = dataset_path('joint-calling/v7/variant_qc/frequencies.ht/', 'analysis')


def filter_lowly_expressed_genes(expression_df):
    """Remove genes with low expression in all samples

    Input:
    expression_df: a data frame with samples as rows and genes as columns, 
    containing normalised expression values (i.e., the average number of molecules 
    for each gene detected in each person).

    Returns:
    A filtered version of the input data frame, after removing columns (genes) 
    with 0 values in more than 10% of the rows (samples).
    """

    # Remove genes with 0 expression in all samples
    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
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


def remove_sc_outliers(df):
    # remove outlier samples, as identified by sc analysis
    outliers = ['966_967', '88_88']
    df = df[-df.sampleid.isin(outliers)]

    return df


def get_number_of_scatters(expression_df, geneloc_df):
    """get index of total number of genes

    Input:
    expression_df: a data frame with samples as rows and genes as columns, 
    containing normalised expression values (i.e., the average number of molecules 
    for each gene detected in each person).
    geneloc_df: a data frame with the gene location (chromosome, start, and 
    end locations as columns) specified for each gene (rows).

    Returns:
    The number of genes (as an int) after filtering for lowly expressed genes.
    This integer number gets fed into the number of scatters to run.
    """

    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = list(expression_df.columns.values)[1:]
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(gene_ids)]
    # subset genes to 10 test genes
    # remove this when running on entire dataset
    genes_to_test = ['ERAP2', 'LAMB1', 'POLR2J2', 'GBP7', 'GSTM3', 'CCDC163P', 'CTLA4', 'AIF1', 'BRCA1', 'IL7', 'HDHD5']
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(genes_to_test)]

    return len(geneloc_df.index)


def get_log_expression(expression_df):

    """get logged expression values

    Input:
    expression_df: a data frame with samples as rows and genes as columns, 
    containing normalised expression values (i.e., the average number of molecules 
    for each gene detected in each person).

    Returns:
    The natural logarithm of expression values, with an offset of 1 to avoid undefined values.
    """

    expression_df = filter_lowly_expressed_genes(expression_df)
    sample_ids = expression_df.iloc[:, 0]
    to_log = expression_df.iloc[:, 1:].columns
    log_expression_df = expression_df[to_log].applymap(lambda x: np.log(x + 1))
    # subset genes to 10 test genes
    # remove this when running on entire dataset
    genes_to_test = ['ERAP2', 'LAMB1', 'POLR2J2', 'GBP7', 'GSTM3', 'CCDC163P', 'CTLA4', 'AIF1', 'BRCA1', 'IL7', 'HDHD5']
    log_expression_df = log_expression_df[log_expression_df.columns.intersection(genes_to_test)]
    log_expression_df.insert(loc=0, column='sampleid', value=sample_ids)
    # remove sc outlier samples
    log_expression_df = remove_sc_outliers(log_expression_df)

    return log_expression_df


def calculate_log_cpm(expression_df, output_prefix):
    """Calculate log cpm for each cell type and chromosome

    Input:
    expression_df: a data frame with samples as rows and genes as columns, 
    containing normalised expression values (i.e., the average number of molecules 
    for each gene detected in each person).

    Returns:
    Counts per million mapped reads
    """

    expression_df = filter_lowly_expressed_genes(expression_df)
    sample_ids = expression_df.iloc[:, 0]
    # remove sampleid column and get log expression
    # this can only be done on integers
    expression_df = expression_df.iloc[:, 1:]
    cpm_df = expression_df.apply(lambda x: (x / sum(x)) * 1000000, axis=0)
    log_cpm = np.log(cpm_df + 1)
    # add sampleids back in
    log_cpm = log_cpm.assign(sampleid=list(sample_ids))
    log_cpm_path = AnyPath(output_prefix) / 'log_cpm.tsv'
    with log_cpm_path.open('w') as fp:
        log_cpm.to_csv(fp, index=False)


def prepare_genotype_info(keys_path, expression_path):
    """Calculate log cpm for each cell type and chromosome

    Input:
    keys_path: path to a tsv file with information on
    OneK1K amd CPG IDs (columns) for each sample (rows).
    expression_path: path to a tsv file of expression values (in transcripts per million)
    with genes as columns and samples as rows.

    Returns:
    Path to a hail matrix table, with rows (alleles) filtered on the following requirements:
    1) biallelic, 2) meets VQSR filters, 3) gene quality score higher than 20,
    4) call rate of 0.8, and 5) variants with MAF < 0.01. Columns (samples) are filtered
    on the basis of having rna-seq expression data, i.e., within the filtered log_expression_df
    """

    init_batch()
    filtered_mt_path = output_path('genotype_table.mt', 'tmp')
    if not hl.hadoop_exists(filtered_mt_path):
        expression_df = pd.read_csv(AnyPath(expression_path), sep='\t')
        log_expression_df = get_log_expression(expression_df)
        mt = hl.read_matrix_table(TOB_WGS)
        mt = hl.experimental.densify(mt)
        # filter to biallelic loci only
        mt = mt.filter_rows(hl.len(mt.alleles) == 2)
        # filter out variants that didn't pass the VQSR filter
        mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
        # VQSR does not filter out low quality genotypes. Filter these out
        mt = mt.filter_entries(mt.GQ <= 20, keep=False)
        # filter out samples with a genotype call rate > 0.8 (as in the gnomAD supplementary paper)
        n_samples = mt.count_cols()
        call_rate = 0.8
        mt = mt.filter_rows(
            hl.agg.sum(hl.is_missing(mt.GT)) > (n_samples * call_rate), keep=False
        )
        # filter out variants with MAF <= 0.01
        ht = hl.read_table(FREQ_TABLE)
        mt = mt.annotate_rows(freq=ht[mt.row_key].freq)
        mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
        # add OneK1K IDs to genotype mt
        sampleid_keys = pd.read_csv(AnyPath(keys_path), sep='\t')
        genotype_samples = pd.DataFrame(mt.s.collect(), columns=['sampleid'])
        sampleid_keys = pd.merge(
            genotype_samples,
            sampleid_keys,
            how='left',
            left_on='sampleid',
            right_on='InternalID',
        )
        sampleid_keys.fillna('', inplace=True)
        sampleid_keys = hl.Table.from_pandas(sampleid_keys)
        sampleid_keys = sampleid_keys.key_by('sampleid')
        mt = mt.annotate_cols(onek1k_id=sampleid_keys[mt.s].OneK1K_ID)
        # only keep samples that have rna-seq expression data
        samples_to_keep = set(log_expression_df.sampleid)
        set_to_keep = hl.literal(samples_to_keep)
        mt = mt.filter_cols(set_to_keep.contains(mt['onek1k_id']))
        mt.write(filtered_mt_path)

    return filtered_mt_path


def calculate_residuals(expression_df, covariate_df, output_prefix):
    """Calculate residuals for each gene in scatter
    
    Input:
    expression_df: a data frame with samples as rows and genes as columns, 
    containing normalised expression values (i.e., the average number of molecules 
    for each gene detected in each person).
    covariate_df: a data frame with samples as rows and cpvariate information 
    as columns (PCA scores, sex, age, and PEER factors).

    Returns: a dataframe of expression residuals, with genes as columns and samples
    as rows.
    """

    log_expression_df = get_log_expression(expression_df)
    gene_ids = list(log_expression_df.columns.values)[1:]
    sample_ids = log_expression_df.merge(covariate_df, on='sampleid', how='right').dropna(axis=0, how='any').sampleid

    # Calculate expression residuals
    def calculate_gene_residual(gene_id):
        """Calculate gene residuals"""
        gene = gene_id
        exprs_val = log_expression_df[['sampleid', gene]]
        test_df = exprs_val.merge(covariate_df, on='sampleid', how='right').dropna(axis=0, how='any')
        test_df = test_df.rename(columns={test_df.columns[1]: 'expression'})
        test_df[['sex', 'age']] = test_df[['sex', 'age']].astype(int)
        y, x = dmatrices(
            'expression ~ sex + PC1 + PC2 + PC3 + PC4 + age + pf1 + pf2',
            test_df,
        )
        model = sm.OLS(y, x)
        residuals = list(model.fit().resid)
        return residuals

    residual_df = pd.DataFrame(list(map(calculate_gene_residual, gene_ids))).T
    residual_df.columns = gene_ids
    residual_df = residual_df.assign(sampleid=list(sample_ids))
    residual_path = AnyPath(output_prefix) / 'log_residuals.tsv'
    with residual_path.open('w') as fp:
        residual_df.to_csv(fp, index=False)

    return residual_df


# Run Spearman rank in parallel by sending genes in batches
def run_spearman_correlation_scatter(
    idx,
    expression,
    geneloc,
    residuals_df,
    filtered_mt_path,
):  # pylint: disable=too-many-locals
    """Run genes in scatter
    
    Input:
    idx: the index of a gene within the filtered expression_df. Note: a 1Mb region is taken 
    up and downstream of each gene.
    expression: input path for the expression data, which consists of a data frame 
    with samples as rows and genes as columns, containing normalised expression 
    values.
    geneloc: input path for gene location data, which consists of a data frame with the 
    gene location (chromosome, start, and end locations as columns) specified for each gene (rows)
    residuals_df: residual dataframe, calculated in calculate_residuals(). 
    filtered_mt_path: the path to a filtered matrix table generated using prepare_genotype_info(), 
    (and filtering requirements outlined in the same function).

    Returns:
    A dataframe with significance results (p-value and FDR) from the Spearman rank correlation. Each
    row contains information for a SNP, and columns contain information on significance levels, 
    spearmans rho, and gene information.
    """

    # calculate log expression values
    expression_df = pd.read_csv(AnyPath(expression), sep='\t')
    log_expression_df = get_log_expression(expression_df)

    # Get 1Mb sliding window around each gene
    geneloc_df = pd.read_csv(AnyPath(geneloc), sep='\t')
    gene_ids = list(log_expression_df.columns.values)[1:]
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(gene_ids)]
    geneloc_df = geneloc_df.assign(left=geneloc_df.start - 1000000)
    geneloc_df = geneloc_df.assign(right=geneloc_df.end + 1000000)
    # subset genes to 10 test genes
    # remove this when running on entire dataset
    genes_to_test = ['ERAP2', 'LAMB1', 'POLR2J2', 'GBP7', 'GSTM3', 'CCDC163P', 'CTLA4', 'AIF1', 'BRCA1', 'IL7', 'HDHD5']
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(genes_to_test)]

    # perform correlation in chunks by gene
    gene_info = geneloc_df.iloc[idx]
    chromosome = gene_info.chr
    # get all SNPs which are within 1Mb of each gene
    init_batch()
    mt = hl.read_matrix_table(filtered_mt_path)
    # only keep samples that are contained within the residuals df
    # this is important, since not all individuals have expression/residual
    # data (this varies by cell type)
    # this also filters out any sc sample outliers
    samples_to_keep = hl.literal(list(residuals_df.sampleid))
    mt = mt.filter_cols(samples_to_keep.contains(mt['onek1k_id']))

    # Do this only on SNPs contained within 1Mb gene region to save on
    # computational time
    first_and_last_snp = chromosome + ':' + str(gene_info.left) + '-' + str(gene_info.right+1)
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(first_and_last_snp, reference_genome='GRCh38')])
    position_table = mt.rows().select()
    position_table = position_table.annotate(
        position=position_table.locus.position,
        snpid=hl.str(position_table.locus)
        + ':'
        + hl.str(position_table.alleles[0])
        + ':'
        + hl.str(position_table.alleles[1]),
    )
    position_table = position_table.to_pandas()
    snps_within_region = position_table[
        position_table['position'].between(gene_info['left'], gene_info['right'])
    ]
    gene_snp_df = snps_within_region.assign(
        gene_id=gene_info.gene_id, gene_symbol=gene_info.gene_name
    )
    
    # get genotypes from mt in order to load individual SNPs into
    # the spearman correlation function
    t = mt.entries()
    t = t.annotate(n_alt_alleles=t.GT.n_alt_alleles())
    t = t.key_by(contig=t.locus.contig, position=t.locus.position)
    t = t.select(t.alleles, t.n_alt_alleles, sampleid=t.onek1k_id)
    t = t.annotate(
        snpid=hl.str(t.contig)
        + ':'
        + hl.str(t.position)
        + ':'
        + hl.str(t.alleles[0])
        + ':'
        + hl.str(t.alleles[1])
    )
    # Do this only on SNPs contained within gene_snp_df to save on
    # computational time
    snps_to_keep = set(gene_snp_df.snpid)
    set_to_keep = hl.literal(snps_to_keep)
    t = t.filter(set_to_keep.contains(t['snpid']))
    # only keep SNPs where all samples have an alt_allele value
    # (argument must be a string, a bytes-like object or a number)
    snps_to_remove = set(t.filter(hl.is_missing(t.n_alt_alleles)).snpid.collect())
    if len(snps_to_remove) > 0:
        t = t.filter(~hl.literal(snps_to_remove).contains(t.snpid))
    
    genotype_df = t.to_pandas(flatten=True)
    # filter gene_snp_df to have the same snps after filtering SNPs that
    # don't have an alt_allele value 
    gene_snp_df = gene_snp_df[gene_snp_df.snpid.isin(set(genotype_df.snpid))]

    # define spearman correlation function, then compute for each SNP
    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_symbol = df.gene_symbol
        gene_id = df.gene_id
        snp = df.snpid
        print(snp)
        gt = genotype_df[genotype_df.snpid == snp][['sampleid', 'n_alt_alleles']]
        res_val = residuals_df[['sampleid', gene_symbol]]
        test_df = res_val.merge(gt, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        # set spearmanr calculation to perform the calculation ignoring nan values
        coef, p = spearmanr(test_df['SNP'], test_df['residual'], nan_policy='omit')
        return (gene_symbol, gene_id, snp, coef, p)

    # run spearman correlation function
    spearman_df = pd.DataFrame(list(gene_snp_df.apply(spearman_correlation, axis=1)))
    print(f'printing spearman df: {spearman_df.head()}')
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
    # turn back into hail table and annotate with global bp,
    # locus, and alleles
    t = hl.Table.from_pandas(spearman_df)
    t = t.annotate(global_bp=hl.locus(t.chrom, hl.int32(t.bp)).global_position())
    t = t.annotate(locus=hl.locus(t.chrom, hl.int32(t.bp)))
    # get alleles
    mt = mt.key_rows_by('locus')
    t = t.key_by('locus')
    t = t.annotate(
        a1=mt.rows()[t.locus].alleles[0],
        a2=hl.if_else(
            hl.len(mt.rows()[t.locus].alleles) == 2, mt.rows()[t.locus].alleles[1], 'NA'
        ),
    )
    # add ID annotation after adding in alleles, a1, and a2
    t = t.annotate(
        id=hl.str(':').join(
            [
                hl.str(t.chrom),
                hl.str(t.bp),
                t.a1,
                t.a2,
                t.gene_symbol,
                # result.db_key, # cell_type_id (eg nk, mononc)
                hl.str(t.round),
            ]
        )
    )
    print(f'Printing final table: {t.show()}')
    spearman_df = t.to_pandas()
    return spearman_df


def merge_df_and_convert_to_string(*df_list):
    """Merge all Spearman dfs and convert to string using .to_string() on df"""
    merged_df: pd.DataFrame = pd.concat(df_list)
    pvalues = merged_df['p_value']
    _, qvals = qvalue(pvalues)
    fdr_values = pd.DataFrame(list(qvals)).iloc[1]
    merged_df = merged_df.assign(fdr=fdr_values)
    merged_df['fdr'] = merged_df.fdr.astype(float)
    return merged_df.to_string()


# Create click command line to enter dependency files
@click.command()
@click.option(
    '--expression', required=True, help='A sample x gene TSV of expression values'
)
@click.option(
    '--geneloc', required=True, help='A TSV of start and end positions for each gene'
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
    geneloc,
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

    # load in files to get number of scatters and residuals
    expression_df_literal = pd.read_csv(AnyPath(expression), sep='\t')
    geneloc_df_literal = pd.read_csv(AnyPath(geneloc), sep='\t')
    covariate_df_literal = pd.read_csv(AnyPath(covariates), sep=',')

    # only run batch jobs if there are test genes in the chromosome
    n_genes_in_scatter = get_number_of_scatters(expression_df_literal, geneloc_df_literal)
    if n_genes_in_scatter == 0:
        raise ValueError('No genes in get_number_of_scatters()')

    calculate_residuals_job = batch.new_python_job('calculate-residuals')
    residuals_df = calculate_residuals_job.call(
        calculate_residuals,
        expression_df=expression_df_literal,
        covariate_df=covariate_df_literal,
        output_prefix=output_prefix,
    )
    # calculate and save log cpm values
    calculate_log_cpm_job = batch.new_python_job('calculate-log-cpm')
    calculate_log_cpm_job.call(
        calculate_log_cpm,
        expression_df=expression_df_literal,
        output_prefix=output_prefix,
    )

    spearman_dfs_from_scatter = []

    filter_mt_job = batch.new_python_job('filter_mt')
    copy_common_env(filter_mt_job)
    filtered_mt_path = filter_mt_job.call(
        prepare_genotype_info, keys_path=keys, expression_path=expression
    )

    for idx in range(n_genes_in_scatter):
        j = batch.new_python_job(name=f'process_{idx}')
        j.cpu(2)
        j.memory('8Gi')
        j.storage('2Gi')
        copy_common_env(j)
        result: hb.resource.PythonResult = j.call(
            run_spearman_correlation_scatter,
            idx=idx,
            expression=expression,
            geneloc=geneloc,
            residuals_df=residuals_df,
            filtered_mt_path=filtered_mt_path,
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
