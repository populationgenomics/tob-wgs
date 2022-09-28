#!/usr/bin/env python3
# flake8: noqa: E402
# pylint: disable=import-outside-toplevel,too-many-locals,import-error,wrong-import-position,too-many-lines,invalid-unary-operand-type,too-many-statements
"""
Create a Hail Batch workflow for all the EQTL analysis, including:

- generate_eqtl_spearman
- conditional_analysis rounds

For example:

    python3 scripts/hail_batch/eqtl_hail_batch/launch_eqtl_spearman.py \
        --input-files-prefix 'gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files' \
        --chromosomes '1 2' \
        --genes B_intermediate
"""
import os
import logging


def setup_logger(name):
    """Get a logger, and silence other common ones"""
    loggers_to_silence = [
        'urllib3.connectionpool',
        'urllib3.util.retry',
        'google.auth._default',
        'google.auth.transport.requests',
        'fsspec',
        'gcsfs',
        'gcsfs.credentials',
    ]
    for logger_name in loggers_to_silence:
        logging.getLogger(logger_name).setLevel(level=logging.WARN)

    logger = logging.getLogger(name)
    logger.setLevel(level=logging.INFO)

    return logger


_logger = setup_logger(os.path.basename(__file__))

import re
from collections import defaultdict
from typing import Any

import click
import hailtop.batch as hb
from cloudpathlib import AnyPath
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    copy_common_env,
    dataset_path,
    init_batch,
    output_path,
    remote_tmpdir,
)
from google.cloud import storage

# computation imports
import hail as hl
import numpy as np
import pandas as pd

DEFAULT_JOINT_CALL_TABLE_PATH = dataset_path('mt/v7.mt/')
DEFAULT_FREQUENCY_TABLE_PATH = dataset_path(
    'joint-calling/v7/variant_qc/frequencies.ht/', 'analysis'
)
DEFAULT_VEP_ANNOTATION_TABLE_PATH = dataset_path('tob_wgs_vep/104/vep104.3_GRCh38.ht/')
DEFAULT_GENCODE_GTF_PATH = 'gs://cpg-reference/gencode/gencode.v38.annotation.gtf.bgz'  # reference_path('gencode_gtf'),

MULTIPY_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16'


# region FILTER_JOINT_CALL_MT


def filter_joint_call_mt(
    *,
    keys_path: str,
    joint_mt_path: str,
    frequency_table_path: str,
    vep_annotation_path: str,
    output_location: str,
    force: bool = False,
):
    """Filter hail matrix table

    Input:
    keys_path: path to a tsv file with information on
    OneK1K amd CPG IDs (columns) for each sample (rows).

    Returns:
    Path to a hail matrix table, with rows (alleles) filtered on the following requirements:
    1) biallelic, 2) meets VQSR filters, 3) gene quality score higher than 20,
    4) call rate of 0.8, and 5) variants with MAF <= 0.01.
    """

    logger = setup_logger('filter_joint_call_mt')
    if AnyPath(output_location).exists() and not force:
        logger.info(f'Reusing existing filtered mt: {output_location}')
        return output_location

    init_batch()
    mt = hl.read_matrix_table(joint_mt_path)
    samples = mt.s.collect()
    n_samples = len(samples)
    mt = mt.naive_coalesce(10000)
    mt = hl.experimental.densify(mt)
    # filter to biallelic loci only
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    # filter out variants that didn't pass the VQSR filter
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    # VQSR does not filter out low quality genotypes. Filter these out
    mt = mt.filter_entries(mt.GQ <= 20, keep=False)
    # filter out samples with a genotype call rate > 0.8 (as in the gnomAD supplementary paper)
    call_rate = 0.8
    mt = mt.filter_rows(
        hl.agg.sum(hl.is_missing(mt.GT)) > (n_samples * call_rate), keep=False
    )
    # filter out variants with MAF <= 0.01
    ht = hl.read_table(frequency_table_path)
    mt = mt.annotate_rows(freq=ht[mt.row_key].freq)
    mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
    # add in VEP annotation
    vep = hl.read_table(vep_annotation_path)
    mt = mt.annotate_rows(
        vep_functional_anno=vep[mt.row_key].vep.regulatory_feature_consequences.biotype
    )
    # add OneK1K IDs to genotype mt
    sampleid_keys = pd.read_csv(AnyPath(keys_path), sep='\t')
    genotype_samples = pd.DataFrame(samples, columns=['sampleid'])
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
    # repartition to save overhead cost
    mt = mt.naive_coalesce(1000)
    mt.write(output_location, overwrite=True)

    return output_location


# endregion FILTER_JOINT_CALL_MT

# region EQTL_SPEARMAN


def generate_eqtl_spearman(
    batch: hb.Batch,
    job_prefix: str,
    *,
    filtered_mt_path: str,
    gencode_gtf_path: str,
    expression_tsv_path: str,
    geneloc_tsv_path: str,
    covariates_tsv_path: str,
    output_prefix: str,
    # can be derived, but better to pass
    genes: list[str],
    gene_level_parallelism: int,
    cell_type: str,
    chromosome: str,  # pylint: disable=unused-argument
    force: bool = False,
):
    """
    Creates a Hail Batch pipeline for calculating EQTLs
    """
    if len(genes) == 0:
        # only run batch jobs if there are test genes in the chromosome
        raise ValueError('No genes for eqtl spearman analysis')

    residuals_csv_path = os.path.join(output_prefix, 'log_residuals.csv')
    if AnyPath(residuals_csv_path).exists() and not force:
        _logger.info(f'Reusing results for eqtl_residuals: {residuals_csv_path}')
    else:
        calculate_residuals_job = batch.new_python_job(
            f'{job_prefix}calculate-residuals'
        )
        copy_common_env(calculate_residuals_job)
        residuals_csv_path = calculate_residuals_job.call(
            calculate_eqtl_residuals,
            output_location=residuals_csv_path,
            expression_tsv_path=expression_tsv_path,
            covariates_tsv_path=covariates_tsv_path,
            force=force,
        )
    # calculate and save log cpm values
    data_summary_path = os.path.join(output_prefix, 'gene_expression.parquet')
    if AnyPath(data_summary_path).exists() and not force:
        _logger.info(f'Reusing results for log_cpm: {data_summary_path}')
    else:
        generate_log_cpm_output_job = batch.new_python_job(
            f'{job_prefix}generate_log_cpm_output'
        )
        generate_log_cpm_output_job.memory('8Gi')
        copy_common_env(generate_log_cpm_output_job)
        generate_log_cpm_output_job.call(
            generate_log_cpm_output,
            cell_type=cell_type,
            output_location=data_summary_path,
            expression_tsv_path=expression_tsv_path,
            gencode_gtf_path=gencode_gtf_path,
            force=force,
        )

    spearman_output_directory = os.path.join(output_prefix, 'parquet')
    jobs = []

    for gene_name in genes:
        spearman_output = os.path.join(
            spearman_output_directory, f'correlation_results_{gene_name}.parquet'
        )

        if AnyPath(spearman_output).exists() and not force:
            _logger.info(
                f'Reusing results for eqtl_spearman_correlation_{gene_name}: {residuals_csv_path}'
            )
            continue

        j = batch.new_python_job(
            name=f'{job_prefix}calculate_spearman_correlation_{gene_name}'
        )
        j.cpu(2)
        j.memory('8Gi')
        j.storage('2Gi')
        j.image(MULTIPY_IMAGE)
        copy_common_env(j)
        j.call(
            run_spearman_correlation_scatter,
            gene_name=gene_name,
            expression_tsv_path=expression_tsv_path,
            geneloc_tsv_path=geneloc_tsv_path,
            residuals_csv_path=residuals_csv_path,
            filtered_mt_path=filtered_mt_path,
            cell_type=cell_type,
            output_location=spearman_output,
            force=force,
        )
        # we'll track job dependencies with an explicit sink
        jobs.append(j)

        # Limit parallelism.
        if len(jobs) >= gene_level_parallelism:
            j.depends_on(jobs[len(jobs) - gene_level_parallelism])

    if len(jobs) == 0:
        # we've reused results the whole way, so no need for 'fake' sink
        sinked_output_dir = spearman_output_directory
    else:
        # create a fake job to avoid having to carry dependencies outside this script
        def return_value(value):
            return value

        fake_sink = batch.new_python_job(f'{job_prefix}sink')
        fake_sink.depends_on(*jobs)
        sinked_output_dir = fake_sink.call(return_value, spearman_output_directory)

    return {
        'spearman_parquet_directory': sinked_output_dir,
        'residuals_csv_path': residuals_csv_path,
    }


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
    """
    Remove outlier samples, as identified by sc analysis
    """
    outliers = ['966_967', '88_88']
    df = df[-df.sampleid.isin(outliers)]

    return df


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
    log_expression_df.insert(loc=0, column='sampleid', value=sample_ids)
    # remove sc outlier samples
    log_expression_df = remove_sc_outliers(log_expression_df)

    return log_expression_df


def calculate_log_cpm(expression_df):
    """
    Remove sampleid column and get log expression
    this can only be done on integers
    """
    expression_df = filter_lowly_expressed_genes(expression_df)
    # remove sampleid column and get log expression
    # this can only be done on integers
    expression_df = expression_df.iloc[:, 1:]
    cpm_df = expression_df.apply(lambda x: (x / sum(x)) * 1000000, axis=0)
    log_cpm = np.log(cpm_df + 1)

    return log_cpm


def generate_log_cpm_output(
    *,
    cell_type: str,
    output_location: str,
    expression_tsv_path: str,
    gencode_gtf_path: str,
    force: bool = False,
):
    """Calculate log cpm for each cell type and chromosome

    Input:
    expression_df: a data frame with samples as rows and genes as columns,
    containing normalised expression values (i.e., the average number of molecules
    for each gene detected in each person).

    Returns:
    Counts per million mapped reads
    """
    logger = setup_logger('generate_log_cpm_output')

    if AnyPath(output_location).exists() and not force:
        logger.info(f'Reusing results from {output_location}')
        return output_location

    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep='\t')
    log_cpm = calculate_log_cpm(expression_df)
    # get gene expression distribution in the output
    # specified here https://github.com/populationgenomics/tob-wgs/issues/97

    def create_struct(gene):
        hist, bin_edges = np.histogram(gene, bins=10)
        n_samples = gene.count()
        min_val = gene.min()
        max_val = gene.max()
        mean_val = gene.mean()
        q1, median_val, q3 = gene.quantile([0.25, 0.5, 0.75])
        iqr = q3 - q1
        # not sure what this value is doing here
        iqr_outliers_min = q1 - (1.5 * iqr)
        iqr_outliers_max = q3 + (1.5 * iqr)
        data_struct = {
            'bin_counts': hist,
            'bin_edges': bin_edges,
            'n_samples': n_samples,
            'min': min_val,
            'max': max_val,
            'mean': mean_val,
            'median': median_val,
            'q1': q1,
            'q3': q3,
            'iqr': iqr,
            'iqr_min': iqr_outliers_min,
            'iqr_max': iqr_outliers_max,
        }

        return data_struct

    # apply function over all columns (genes) and reset index so that gene names are in the df
    data_summary = pd.DataFrame(
        log_cpm.apply(create_struct, axis=0), columns=['data']
    ).reset_index()
    data_summary = data_summary.rename({'index': 'gene_name'}, axis='columns')
    # add in cell type info
    data_summary['cell_type_id'] = cell_type
    # add in ENSEMBL IDs
    init_batch(driver_cores=8)
    gtf = hl.experimental.import_gtf(
        gencode_gtf_path, reference_genome='GRCh38', skip_invalid_contigs=True
    )
    # convert int to str in order to avoid 'int() argument must be a string, a bytes-like object or a number, not 'NoneType''
    gtf = gtf.annotate(frame=hl.str(gtf.frame))
    gtf = gtf.to_pandas()
    data_summary['ensembl_ids'] = data_summary.merge(
        gtf.drop_duplicates('gene_name'), how='left', on='gene_name'
    ).gene_id

    # Save file
    with AnyPath(output_location).open('wb+', force_overwrite_to_cloud=True) as fp:
        data_summary.to_parquet(fp)

    return output_location


def calculate_eqtl_residuals(
    expression_tsv_path: str,
    covariates_tsv_path: str,
    output_location: str,
    force: bool = False,
):
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

    logger = setup_logger('calculate_residuals')

    if AnyPath(output_location).exists() and not force:
        logger.info(f'Reusing previous results from {output_location}')
        return output_location

    # import these here to avoid polluting the global namespace
    #   (and make testing easier)
    import statsmodels.api as sm
    from patsy import dmatrices  # pylint: disable=no-name-in-module

    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep='\t')
    covariates_df = pd.read_csv(AnyPath(covariates_tsv_path), sep=',')

    log_expression_df = get_log_expression(expression_df)
    gene_ids = list(log_expression_df.columns.values)[1:]
    sample_ids = (
        log_expression_df.merge(covariates_df, on='sampleid', how='right')
        .dropna(axis=0, how='any')
        .sampleid
    )

    # Calculate expression residuals
    def calculate_gene_residual(gene_id):
        """Calculate gene residuals"""
        gene = gene_id
        exprs_val = log_expression_df[['sampleid', gene]]
        test_df = exprs_val.merge(covariates_df, on='sampleid', how='right').dropna(
            axis=0, how='any'
        )
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

    with AnyPath(output_location).open('w+', force_overwrite_to_cloud=True) as fp:
        residual_df.to_csv(fp, index=False)

    return output_location


# Run Spearman rank in parallel by sending genes in batches
def run_spearman_correlation_scatter(
    *,
    expression_tsv_path,
    geneloc_tsv_path,
    residuals_csv_path,
    filtered_mt_path,
    cell_type,
    gene_name: str,
    output_location: str,
    force: bool = False,
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

    # checkpoint
    if AnyPath(output_location).exists() and not force:
        return output_location

    # silences GCSFS messages
    _ = setup_logger('run_spearman_correlation_scatter')

    # import multipy here to avoid issues with driver image updates
    from multipy.fdr import qvalue
    from scipy.stats import spearmanr

    # calculate log expression values
    residuals_df = pd.read_csv(residuals_csv_path)
    expression_df = pd.read_csv(expression_tsv_path, sep='\t')
    # log_expression_df = get_log_expression(expression_df)

    # Get 1Mb sliding window around each gene
    geneloc_df = pd.read_csv(geneloc_tsv_path, sep='\t')
    gene_infos = geneloc_df[geneloc_df.gene_name == gene_name]
    if len(gene_infos) == 0:
        raise ValueError(
            f'Could not find gene {gene_name} in geneloc_df: {geneloc_tsv_path}'
        )
    gene_info = gene_infos.iloc[0]
    gene_info.left = int(gene_info.start - 1e6)
    gene_info.right = int(gene_info.end + 1e6)
    chromosome = gene_info.chr

    # perform correlation in chunks by gene

    # get all SNPs which are within 1Mb of each gene
    # TODO: can we drop these cores: 8 -> 1
    init_batch(driver_cores=1, worker_cores=1, driver_memory='highmem')
    mt = hl.read_matrix_table(filtered_mt_path)
    # only keep samples that are contained within the residuals df
    # this is important, since not all individuals have expression/residual
    # data (this varies by cell type)
    # this also filters out any sc sample outliers
    samples_to_keep = hl.literal(list(residuals_df.sampleid))
    # Do this only on SNPs contained within 1Mb gene region to save on
    # computational time
    right_boundary = min(
        gene_info.right + 1, hl.get_reference('GRCh38').lengths[chromosome]
    )
    left_boundary = max(1, gene_info.left)
    first_and_last_snp = f'{chromosome}:{left_boundary}-{right_boundary}'
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(first_and_last_snp, reference_genome='GRCh38')]
    )
    # remove SNPs with no variance
    mt = mt.filter_rows(
        (hl.agg.all(mt.GT.n_alt_alleles() == 0))
        | (hl.agg.all(mt.GT.n_alt_alleles() == 1))
        | (hl.agg.all(mt.GT.n_alt_alleles() == 2)),
        keep=False,
    )
    mt = mt.filter_cols(samples_to_keep.contains(mt['onek1k_id']))
    # TODO: mfranklin to check mt.persist, as mt.rows() is called twice
    mt = mt.persist('MEMORY_AND_DISK')

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
        position_table['position'].between(gene_info.left, gene_info.right)
    ]
    gene_snp_df = snps_within_region.assign(
        gene_id=gene_info.gene_id, gene_symbol=gene_info.gene_name
    )

    # get genotypes from mt in order to load individual SNPs into
    # the spearman correlation function
    t = mt.entries()
    t = t.key_by()
    t = t.select(
        n_alt_alleles=t.GT.n_alt_alleles(),
        sampleid=t.onek1k_id,
        contig=t.locus.contig,
        position=t.locus.position,
        global_bp=t.locus.global_position(),
        a1=t.alleles[0],
        a2=t.alleles[1],
    )
    t = t.annotate(
        snpid=hl.str(t.contig)
        + ':'
        + hl.str(t.position)
        + ':'
        + hl.str(t.a1)
        + ':'
        + hl.str(t.a2)
    )
    # Do this only on SNPs contained within gene_snp_df to save on
    # computational time
    snps_to_keep = hl.literal(set(gene_snp_df.snpid))
    t = t.filter(snps_to_keep.contains(t['snpid']))
    # only keep SNPs where all samples have an alt_allele value
    # (argument must be a string, a bytes-like object or a number)
    snps_to_remove = set(t.filter(hl.is_missing(t.n_alt_alleles)).snpid.collect())
    if len(snps_to_remove) > 0:
        t = t.filter(~hl.literal(snps_to_remove).contains(t.snpid))

    genotype_df = t.to_pandas(flatten=True)
    # filter gene_snp_df to have the same snps after filtering SNPs that
    # don't have an alt_allele value
    gene_snp_df = gene_snp_df[gene_snp_df.snpid.isin(set(genotype_df.snpid))]

    # Get association effect data, to be used for violin plots for each genotype of each SNP
    # TODO: remove redefined gene function param
    def get_association_effect_data(gene):
        sampleid = expression_df.sampleid
        log_cpm = calculate_log_cpm(expression_df)
        log_cpm['sampleid'] = sampleid
        # reduce log_cpm matrix to gene of interest only
        log_cpm = log_cpm[['sampleid', gene]]
        # merge log_cpm data with genotype info
        gt_expr_data = genotype_df.merge(
            log_cpm, left_on='sampleid', right_on='sampleid'
        )
        # group gt_expr_data by snpid and genotype
        grouped_gt = gt_expr_data.groupby(['snpid', 'n_alt_alleles'])

        def create_struct(gene, group):
            hist, bin_edges = np.histogram(group[gene], bins=10)
            n_samples = group[gene].count()
            min_val = group[gene].min()
            max_val = group[gene].max()
            mean_val = group[gene].mean()
            q1, median_val, q3 = group[gene].quantile([0.25, 0.5, 0.75])
            iqr = q3 - q1
            data_struct = {
                'bin_counts': hist,
                'bin_edges': bin_edges,
                'n_samples': n_samples,
                'min': min_val,
                'max': max_val,
                'mean': mean_val,
                'median': median_val,
                'q1': q1,
                'q3': q3,
                'iqr': iqr,
            }

            return data_struct

        snp_gt_summary_data = []
        for item in grouped_gt.groups:
            snp_id = item[0]
            genotype = item[1]
            group = grouped_gt.get_group((snp_id, genotype))
            chrom = group.contig.iloc[0]
            bp = group.position.iloc[0]
            global_bp = group.global_bp.iloc[0]
            a1 = group.a1.iloc[0]
            a2 = group.a2.iloc[0]
            gene_id = gene_info.gene_id
            snp_gt_summary_data.append(
                {
                    'chrom': chrom,
                    'bp': bp,
                    'a1': a1,
                    'a2': a2,
                    'global_bp': global_bp,
                    'genotype': genotype,
                    'gene_id': gene_id,
                    'gene_symbol': gene,
                    'cell_type_id': cell_type,
                    'struct': create_struct(gene, group),
                }
            )

        return pd.DataFrame.from_dict(snp_gt_summary_data)

    association_effect_data = get_association_effect_data(gene_name)
    # write association effect data
    association_effect_path = output_path(
        os.path.join(
            'association-effect-data',
            cell_type,
            str(chromosome),
            f'eqtl_effect_{gene_name}.parquet',
        )
    )
    with AnyPath(association_effect_path).open(
        'wb', force_overwrite_to_cloud=True
    ) as fp:
        association_effect_data.to_parquet(fp)

    # define spearman correlation function, then compute for each SNP
    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_symbol = df.gene_symbol
        gene_id = df.gene_id
        alleles = df.alleles
        snp = df.snpid
        gt = genotype_df[genotype_df.snpid == snp][['sampleid', 'n_alt_alleles']]
        res_val = residuals_df[['sampleid', gene_symbol]]
        test_df = res_val.merge(gt, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        # set spearmanr calculation to perform the calculation ignoring nan values
        coef, p = spearmanr(test_df['SNP'], test_df['residual'], nan_policy='omit')
        return gene_symbol, gene_id, alleles, snp, coef, p

    # run spearman correlation function
    spearman_df = pd.DataFrame(list(gene_snp_df.apply(spearman_correlation, axis=1)))
    spearman_df.columns = [
        'gene_symbol',
        'gene_id',
        'alleles',
        'snp_id',
        'spearmans_rho',
        'p_value',
    ]

    # add in locus and chromosome information to get global position in hail
    locus = spearman_df.snp_id.str.split(':', expand=True)[[0, 1]].agg(':'.join, axis=1)
    chrom = locus.str.split(':', expand=True)[0]
    bp = locus.str.split(':', expand=True)[1]
    spearman_df['locus'], spearman_df['chrom'], spearman_df['bp'] = [
        locus,
        chrom,
        bp,
    ]
    # turn into hail table and annotate with global bp and allele info
    t = hl.Table.from_pandas(spearman_df)
    t = t.annotate(global_bp=hl.locus(t.chrom, hl.int32(t.bp)).global_position())
    t = t.annotate(locus=hl.locus(t.chrom, hl.int32(t.bp)))
    t = t.annotate(a1=t.alleles[0], a2=t.alleles[1])
    # add in vep annotation
    t = t.key_by(t.locus, t.alleles)
    t = t.annotate(functional_annotation=mt.rows()[t.key].vep_functional_anno)
    t = t.key_by()
    t = t.drop(t.locus, t.alleles, t.snp_id)
    # turn back into pandas df and add additional information
    # for front-end analysis
    spearman_df = t.to_pandas()
    spearman_df['round'] = '1'
    # add celltype id
    celltype_id = cell_type.lower()
    spearman_df['cell_type_id'] = celltype_id
    # Correct for multiple testing using Storey qvalues
    # qvalues are used instead of BH/other correction methods, as they do not assume independence (e.g., high LD)
    pvalues = spearman_df['p_value']
    _, qvals = qvalue(pvalues)
    fdr_values = pd.DataFrame(list(qvals))
    spearman_df = spearman_df.assign(fdr=fdr_values)
    spearman_df['fdr'] = spearman_df.fdr.astype(float)

    # Save file
    with AnyPath(output_location).open('wb+', force_overwrite_to_cloud=True) as fp:
        spearman_df.to_parquet(fp)

    return output_location


# endregion EQTL_SPEARMAN

# region CONDITIONAL_ANALYSIS


def generate_conditional_analysis(
    batch: hb.Batch,
    *,
    job_prefix: str,
    cell_type: str,
    chromosome: str,  # pylint: disable=unused-argument
    filtered_matrix_table_path: str,
    significant_snps_directory: str,
    residuals_path: str,
    genes: list[str],
    gene_level_parallelism: int,
    iterations: int = 4,
    output_prefix: str,
    force: bool = False,
):
    """
    Creates a Hail Batch pipeline for calculating eQTLs using {iterations} iterations,
    scattered across the number of genes. Note, iteration 1 is completed in `generate_eqtl_spearan.py`.
    """

    # these are now PATHS, that we'll load and unload each bit
    previous_sig_snps_directory = (
        significant_snps_directory  # pylint: disable=invalid-name
    )
    previous_residual_path = residuals_path  # pylint: disable=invalid-name

    sig_snps_dependencies: list[Any] = []

    # Perform conditional analysis for n iterations (specified in click interface)
    # note, iteration 1 is performed in generate_eqtl_spearman.py, which requires starting
    # the iteration at 2 (from a 0 index)
    for iteration in range(2, iterations + 2):
        round_dir = os.path.join(output_prefix, f'round{iteration}/')
        residuals_output_path = os.path.join(round_dir, f'residual_results.csv')

        if AnyPath(residuals_output_path).exists() and not force:
            _logger.info(f'Reusing residuals result from {residuals_output_path}')
            previous_residual_path = residuals_output_path
        else:
            calc_resid_df_job = batch.new_python_job(
                f'{job_prefix}calculate-residual-round-{iteration}'
            )
            calc_resid_df_job.depends_on(*sig_snps_dependencies)
            calc_resid_df_job.cpu(2)
            calc_resid_df_job.memory('8Gi')
            calc_resid_df_job.storage('2Gi')
            copy_common_env(calc_resid_df_job)

            # this calculate_residual_df will also write an output that is useful later
            # if we assign previous_residual_path to the result of the call, we get the job dependency we actually want
            # as opposed to saying the following, which means there's actually no job dependency (bad):
            #    previous_residual_path = <output_path>
            previous_residual_path = calc_resid_df_job.call(
                calculate_conditional_residuals,
                residual_path=previous_residual_path,
                significant_snps_path=previous_sig_snps_directory,
                output_location=os.path.join(round_dir, f'residual_results.csv'),
                filtered_matrix_table_path=filtered_matrix_table_path,
                force=force,
            )

        new_sig_snps_directory = os.path.join(round_dir, 'sigsnps/')
        jobs = []
        for gene_name in genes:
            conditional_output_path = os.path.join(
                new_sig_snps_directory, f'sig-snps-{gene_name}.parquet'
            )

            if AnyPath(conditional_output_path).exists() and not force:
                _logger.info(
                    f'Reusing conditionals result from {residuals_output_path}'
                )
            else:
                j = batch.new_python_job(
                    name=f'{job_prefix}calculate_spearman_round_{iteration}_{gene_name}'
                )
                j.cpu(2)
                j.memory('8Gi')
                j.storage('2Gi')
                j.image(MULTIPY_IMAGE)
                copy_common_env(j)
                j.call(
                    run_scattered_conditional_analysis,
                    iteration=iteration,
                    cell_type=cell_type,
                    gene_name=gene_name,
                    residual_path=previous_residual_path,
                    significant_snps_path=previous_sig_snps_directory,
                    filtered_matrix_table_path=filtered_matrix_table_path,
                    round_outputdir=round_dir,
                    output_location=conditional_output_path,
                    force=force,
                )
                if sig_snps_dependencies:
                    j.depends_on(*sig_snps_dependencies)
                jobs.append(j)

                # Limit parallelism.
                if len(jobs) >= gene_level_parallelism:
                    j.depends_on(jobs[len(jobs) - gene_level_parallelism])

        previous_sig_snps_directory = new_sig_snps_directory
        sig_snps_dependencies = jobs

    if len(sig_snps_dependencies) == 0:
        # we've reused results the whole way, so no need for 'fake' sink
        sinked_sig_snps_directory = previous_sig_snps_directory
    else:
        # create a fake job to avoid having to carry dependencies outside this script
        def return_value(value):
            return value

        fake_sink = batch.new_python_job(f'{job_prefix}sink')
        fake_sink.depends_on(*sig_snps_dependencies)
        sinked_sig_snps_directory = fake_sink.call(
            return_value, previous_sig_snps_directory
        )

    return {
        'sig_snps_parquet_directory': sinked_sig_snps_directory,
        'residuals_csv_path': previous_residual_path,
    }


def get_genotype_df(residual_df, gene_snp_test_df, filtered_matrix_table_path):
    """load genotype df and filter
    Input:
    residual_df: residual dataframe, calculated in calculate_residual_df(). This is run for
    each round/iteration of conditional analysis. The dataframe consists of genes as columns
    and samples as rows.
    gene_snp_test_df: a dataframe with rows as SNPs and columns containing snp_id, gene_symbol, and
    gene_id information.

    Returns:
    A pandas dataframe, with genotypes for each sample (in rows) for every SNP. Genotype information is
    stored as the number of alternative alleles (n_alt_alleles; 0, 1, or 2).
    """

    init_batch()
    mt = hl.read_matrix_table(filtered_matrix_table_path)
    # only keep samples that are contained within the residuals df
    # this is important, since not all indivuduals have expression/residual
    # data (this varies by cell type)
    samples_to_keep = hl.set(list(residual_df.sampleid))

    # Do this only on SNPs contained within gene_snp_df to save on
    # computational time
    snps_to_keep = set(gene_snp_test_df.snp_id)
    sorted_snps = sorted(snps_to_keep)
    sorted_snp_positions = list(map(lambda x: x.split(':')[:2][1], sorted_snps))
    sorted_snp_positions = [int(i) for i in sorted_snp_positions]
    # get first and last positions, with 1 added to last position (to make it inclusive)
    # get one random element from the list to get chromosome
    random_snp = next(iter(snps_to_keep))
    chromosome = random_snp.split(':')[0]

    first_and_last_snp = (
        chromosome
        + ':'
        + str(sorted_snp_positions[0])
        + '-'
        + str(sorted_snp_positions[-1] + 1)
    )
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(first_and_last_snp, reference_genome='GRCh38')]
    )
    mt = mt.filter_cols(samples_to_keep.contains(mt['onek1k_id']))
    t = mt.entries()
    t = t.annotate(n_alt_alleles=t.GT.n_alt_alleles())
    t = t.key_by(contig=t.locus.contig, position=t.locus.position)
    t = t.select(t.alleles, t.n_alt_alleles, sampleid=t.onek1k_id)
    t = t.annotate(
        snp_id=hl.str(':').join(
            list(map(hl.str, [t.contig, t.position, t.alleles[0], t.alleles[1]]))
        )
    )
    # Further reduce the table by only selecting SNPs needed
    set_to_keep = hl.literal(snps_to_keep)
    t = t.filter(set_to_keep.contains(t['snp_id']))
    # only keep SNPs where all samples have an alt_allele value
    snps_to_remove = set(t.filter(hl.is_missing(t.n_alt_alleles)).snp_id.collect())
    if len(snps_to_remove) > 0:
        t = t.filter(~hl.literal(snps_to_remove).contains(t.snp_id))

    genotype_df = t.to_pandas(flatten=True)
    genotype_df.rename({'onek1k_id': 'sampleid'}, axis=1, inplace=True)

    return genotype_df


def calculate_conditional_residuals(
    residual_path: str,
    significant_snps_path: str,
    output_location: str,
    filtered_matrix_table_path: str,
    force: bool = False,
) -> str:
    """calculate residuals for gene list

    Input:
    residual_df: a dataframe of expression residuals, with genes as columns and samples
    as rows.
    significant_snps_df: a dataframe with significance results (p-value and FDR) from the Spearman rank
    correlation calculated in `conditional_analysis.py`. Each row contains information for a SNP,
    and columns contain information on significance levels, spearmans rho, and gene information.

    Returns: a dataframe of expression residuals, with genes as columns and samples
    as rows, which have been conditioned on the lead eSNP.
    """

    # set default logger to ERROR only

    logger = setup_logger('calculate_residual_df')

    logger.info(f'Got paths: {residual_path} / {significant_snps_path}')

    if AnyPath(output_location).exists() and not force:
        logger.info(f'Reusing residuals result from {output_location}')
        return output_location

    # import these here to avoid polluting the global namespace
    #   (and make testing easier)
    import statsmodels.api as sm
    from patsy import dmatrices  # pylint: disable=no-name-in-module

    residual_df = pd.read_csv(residual_path)
    significant_snps_df = pd.read_parquet(significant_snps_path)

    # make sure 'gene_symbol' is the first column
    # otherwise, error thrown when using reset_index
    cols = list(significant_snps_df)
    cols.insert(0, cols.pop(cols.index('gene_symbol')))
    significant_snps_df = significant_snps_df.loc[:, cols]

    # Identify the top eSNP for each eGene
    esnp1 = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .first()
        .reset_index()
    )
    # add in SNP_ID column
    esnp1['snp_id'] = esnp1[['chrom', 'bp', 'a1', 'a2']].apply(
        lambda row: ':'.join(row.values.astype(str)), axis=1
    )
    # Subset residuals for the genes to be tested
    gene_ids = esnp1['gene_symbol'][
        esnp1['gene_symbol'].isin(residual_df.columns)
    ].to_list()
    # save sampleids before filtering redidual_df
    sample_ids = residual_df.sampleid
    residual_df = residual_df[gene_ids]
    # reassign sample ids
    residual_df = residual_df.assign(sampleid=sample_ids.to_list())
    # Subset genotype file for the significant SNPs
    genotype_df = get_genotype_df(
        residual_df=residual_df,
        gene_snp_test_df=esnp1,
        filtered_matrix_table_path=filtered_matrix_table_path,
    )

    # Find residuals after adjustment of lead SNP
    def calculate_adjusted_residuals(gene_id):
        # select gene to regress
        exprs_val = residual_df[['sampleid', gene_id]]
        # In selecting the top esnp, let's double check the gene is in the esnp1
        esnps_for_gene = esnp1.snp_id[esnp1.gene_symbol == gene_id].values
        if len(esnps_for_gene) == 0:
            # throw an error if the gene (for some reason) isn't in esnp1
            raise ValueError(f'Selected gene {gene_id} was not found in esnp1')
        # get top esnp for the gene
        snp = esnps_for_gene[0]
        snp_genotype = genotype_df[genotype_df.snp_id == snp][
            ['sampleid', 'n_alt_alleles']
        ]

        # Create a test df by adding covariates
        test_df = exprs_val.merge(snp_genotype, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'expression', 'genotype']

        y, x = dmatrices('expression ~ genotype', test_df)
        try:
            model = sm.OLS(y, x)
            residuals = list(model.fit().resid)
            return residuals
        except Exception as e:
            logger.info(f'y = {y}, x = {x}')
            raise Exception(
                f'Error during calculate_adjusted_residuals for {gene_id}'
            ) from e

    adjusted_residual_mat = pd.DataFrame(
        list(map(calculate_adjusted_residuals, gene_ids))
    ).T
    adjusted_residual_mat.columns = gene_ids
    adjusted_residual_mat.insert(loc=0, column='sampleid', value=genotype_df.sampleid)

    adjusted_residual_mat.to_csv(output_location, index=False)
    return output_location


def run_scattered_conditional_analysis(
    iteration: int,
    cell_type: str,
    gene_name: int,
    residual_path: str,
    significant_snps_path: str,
    filtered_matrix_table_path: str,
    round_outputdir: str,
    output_location: str,
    force: bool = False,
):
    """Run genes in scatter

    Input:
    iteration: synonymous with round. 4 iterations are performed by default for the conditional
    analysis, which results in 5 total rounds with the addition of the first analysis.
    idx: the index of each gene with a top eSNP. Note: a 1Mb region is taken
    up and downstream of each gene.
    residual_df: residual dataframe, calculated in calculate_residual_df(). This is run for
    each round/iteration of conditional analysis.

    Returns:
    A dataframe with significance results (p-value and FDR) from the Spearman rank correlation. Each
    row contains information for a SNP, and columns contain information on significance levels,
    spearmans rho, and gene information.
    """

    logger = setup_logger('run_scattered_conditional_analysis')

    if AnyPath(output_location).exists() and not force:
        logger.info(f'Reuse conditional output {output_location}')
        return output_location

    # import these here to avoid polluting the global namespace
    from multipy.fdr import qvalue
    from scipy.stats import spearmanr

    residual_df = pd.read_csv(residual_path)
    significant_snps_df = pd.read_parquet(significant_snps_path)
    significant_snps_df = significant_snps_df[
        significant_snps_df.gene_symbol == gene_name
    ]

    if len(significant_snps_df) == 0:
        # previous round failed, because this gene is NOT in the previous significant_snps
        logger.warning('This is skipped, because there is no data from previous rounds')
        return None

    # make sure 'gene_symbol' is the first column
    # otherwise, error thrown when using reset_index
    cols = list(significant_snps_df)
    cols.insert(0, cols.pop(cols.index('gene_symbol')))
    significant_snps_df = significant_snps_df.loc[:, cols]

    # Identify the top eSNP for each eGene
    esnp1 = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .first()
        .reset_index()
    )
    # save esnp1 for front-end use on which SNPs have been conditioned on
    esnp1_path = (
        AnyPath(round_outputdir) / f'conditioned_esnps_{iteration}_{gene_name}.tsv'
    )
    # if this succeeds, and the following code fails, then this will fail with
    # an 'cloudpathlib.exceptions.OverwriteNewerCloudError', so force overwrite
    with esnp1_path.open('w+', force_overwrite_to_cloud=True) as fp:
        esnp1.to_csv(fp, index=False)

    # Remaning eSNPs to test
    esnps_to_test = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .apply(lambda group: group.iloc[1:, 1:])
        .reset_index()
    )
    # add in SNP_ID column
    esnps_to_test = esnps_to_test[esnps_to_test.gene_symbol.isin(residual_df.columns)]
    esnps_to_test['snp_id'] = esnps_to_test[['chrom', 'bp', 'a1', 'a2']].apply(
        lambda row: ':'.join(row.values.astype(str)), axis=1
    )

    # for each gene, get esnps_to_test
    gene_snp_test_df = esnps_to_test[['snp_id', 'gene_symbol', 'gene_id', 'a1', 'a2']]
    gene_snp_test_df = gene_snp_test_df[gene_snp_test_df['gene_symbol'] == gene_name]

    # Subset genotype file for the significant SNPs
    genotype_df = get_genotype_df(
        residual_df=residual_df,
        gene_snp_test_df=gene_snp_test_df,
        filtered_matrix_table_path=filtered_matrix_table_path,
    )
    res_val = residual_df[['sampleid', gene_name]]

    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_id = df.gene_id
        a1 = df.a1
        a2 = df.a2
        snp = df.snp_id
        gt = genotype_df[genotype_df.snp_id == snp][['sampleid', 'n_alt_alleles']]

        test_df = res_val.merge(gt, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        # Remove SNPs with NA values. This is due to some individuals having NA genotype values when
        # the residual df was created using the top eSNP as the conditioning SNP.
        test_df = test_df.dropna(axis=0, how='any')
        spearmans_rho, p = spearmanr(
            test_df['SNP'], test_df['residual'], nan_policy='omit'
        )
        return gene_name, gene_id, a1, a2, snp, spearmans_rho, p

    # calculate spearman correlation
    adjusted_spearman_df = pd.DataFrame(
        list(gene_snp_test_df.apply(spearman_correlation, axis=1))
    )

    adjusted_spearman_df.columns = [
        'gene_symbol',
        'gene_id',
        'a1',
        'a2',
        'snp_id',
        'spearmans_rho',
        'p_value',
    ]
    # remove any NA values (i.e., individuals with zero variance in their genotypes)
    adjusted_spearman_df = adjusted_spearman_df.dropna(axis=0, how='any')
    if len(adjusted_spearman_df) == 0:
        # This is an error case, the correlation coefficient is not defined
        # for any SNPs in this gene. So instead, lets:
        #   - Not write anything to the output directory (filtered out for next round)
        #   - Write a file somewhere to easily track that this failed round
        logger.error(f'All residuals for {gene_name} are 0, returning nothing')
        failed_marker = output_path(
            os.path.join('failed', os.path.basename(output_location))
        )
        with AnyPath(failed_marker).open(mode='w+', force_overwrite_to_cloud=True) as f:
            f.write(f'All residuals for {gene_name} are 0')
        return None

    # add in locus and chromosome information to get global position in hail
    locus = adjusted_spearman_df.snp_id.str.split(':', expand=True)[[0, 1]].agg(
        ':'.join, axis=1
    )
    chrom = locus.str.split(':', expand=True)[0]
    bp = locus.str.split(':', expand=True)[1]
    (
        adjusted_spearman_df['locus'],
        adjusted_spearman_df['chrom'],
        adjusted_spearman_df['bp'],
    ) = [
        locus,
        chrom,
        bp,
    ]
    adjusted_spearman_df['round'] = str(iteration)
    # turn into hail table and annotate with global bp and allele info
    t = hl.Table.from_pandas(adjusted_spearman_df)
    t = t.annotate(global_bp=hl.locus(t.chrom, hl.int32(t.bp)).global_position())
    # turn back into pandas df and add additional information
    # for front-end analysis
    adjusted_spearman_df = t.to_pandas()
    # add celltype id
    celltype_id = cell_type.lower()
    adjusted_spearman_df['cell_type_id'] = celltype_id
    # Correct for multiple testing using Storey qvalues
    # qvalues are used instead of BH/other correction methods, as they do not assume independence (e.g., high LD)
    pvalues = adjusted_spearman_df['p_value']
    _, qvals = qvalue(pvalues)
    fdr_values = pd.DataFrame(list(qvals))
    adjusted_spearman_df = adjusted_spearman_df.assign(fdr=fdr_values)
    adjusted_spearman_df['fdr'] = adjusted_spearman_df.fdr.astype(float)
    adjusted_spearman_df['cell_type_id'] = cell_type

    # Save file
    adjusted_spearman_df.to_parquet(output_location)

    return output_location


# Run Spearman rank in parallel by sending genes in batches


# endregion CONDITIONAL_ANALYSIS

# region GENERIC


def list_dir(directory: str, filter_=None):
    """Generic list_dir implementation with ability to filter"""
    if directory.startswith('gs://'):
        cs_client = storage.Client()
        bucket_name, bucket_path = directory.split('gs://')[1].split('/', maxsplit=1)
        bucket = cs_client.get_bucket(bucket_name)

        blobs = bucket.list_blobs(prefix=bucket_path + '/', delimiter='/')
        return [b.name for b in blobs if not filter_ or filter_(b.name)]

    if os.path.exists(directory):
        return [fn for fn in os.listdir(directory) if not filter_ or filter_(fn)]

    raise ValueError(f'Unrecognised directory type: {directory}')


def get_chromosomes(input_files_prefix: str):
    """
    Get all chromosomes for which gene_location_files exist
    """
    pattern = re.compile(r'GRCh38_geneloc_chr(.+)\.tsv$')
    files = list_dir(
        os.path.join(input_files_prefix, 'gene_location_files'),
        # check if the pattern matches
        filter_=pattern.search,
    )
    searches = [pattern.search(fn) for fn in files]
    chromosomes = set(
        file_search.groups()[0] for file_search in searches if file_search
    )
    return sorted(chromosomes)


def get_cell_types_from(input_files_prefix: str) -> list[str]:
    """
    we can infer the cell types from the 'covariates_files' subdirectory
        eg: {cell_type}_peer_factors_file.tsv
    """
    covariates_files_dir = os.path.join(input_files_prefix, 'covariates_files')
    _logger.info(f'Going to fetch cell types from {covariates_files_dir}')
    ending = '_peer_factors_file.txt'
    return [
        os.path.basename(fn)[: -len(ending)]
        for fn in list_dir(covariates_files_dir, lambda el: el.endswith(ending))
    ]


def get_genes_for_chromosome(*, expression_tsv_path, geneloc_tsv_path) -> list[str]:
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
    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep='\t')
    geneloc_df = pd.read_csv(AnyPath(geneloc_tsv_path), sep='\t')

    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = set(list(expression_df.columns.values)[1:])

    genes = set(geneloc_df.gene_name).intersection(gene_ids)
    return list(sorted(genes))


# endregion GENERIC

# region DRIVER


@click.command()
@click.option(
    '--input-files-prefix',
    required=True,
    help='A path prefix of where input files are located. eg: gs://MyBucket/folder. '
    'If a relative path is given, it will be from the output-path',
)
@click.option(
    '--chromosomes',
    help='List of chromosome numbers to run eQTL analysis on. '
    'Space separated, as one argument (Default: all)',
)
@click.option(
    '--cell-types',
    default=None,
    multiple=True,
    help='List of cell types to test. All available cell types can be found in '
    '`gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/expression_files/`',
)
@click.option(
    '--gene',
    type=str,
    multiple=True,
    help='Limit to specific genes to test with (name of gene)',
)
@click.option(
    '--gene-level-parallelism',
    type=int,
    default=10,
    help='The number of genes to process in parallel, for each chromosome and cell type. '
    'This is useful to avoid starting too many jobs near the root of the scheduling tree, '
    'leaving too few resources for leaf nodes (particularly Hail Query jobs). The default '
    'value is fairly small, as we also process chromosomes and cell types in parallel.',
)
@click.option('--force', is_flag=True, help='Skip checkpoints')
@click.option('--local-debug', is_flag=True, help='Dry run without service-backend')
def from_cli(
    chromosomes: str,
    input_files_prefix: str,
    cell_types: list[str] | None,
    force: bool = False,
    local_debug: bool = False,
    gene: list[str] = None,
    gene_level_parallelism,
):
    """Run the EQTL analysis from command line arguments"""
    chromosomes_list = chromosomes.split(' ') if chromosomes else None
    if local_debug:
        backend = None
        run_args = {'dry_run': True}
    else:
        backend = hb.ServiceBackend(
            billing_project=get_config()['hail']['billing_project'],
            remote_tmpdir=remote_tmpdir(),
        )
        run_args = {'wait': False}

    batch = hb.Batch(
        name='eqtl_spearman', backend=backend, default_python_image=MULTIPY_IMAGE
    )
    main(
        batch,
        input_files_prefix=input_files_prefix,
        chromosomes=chromosomes_list,
        cell_types=cell_types,
        genes=gene,
        gene_level_parallelism=gene_level_parallelism,
        force=force,
    )
    # pylint: disable=protected-access
    _logger.info(f'Got {len(batch._jobs)} jobs in {batch.name}')
    batch.run(**run_args)


def main(
    batch: hb.Batch,
    *,
    input_files_prefix: str,
    chromosomes: list[str],
    cell_types: list[str],
    genes: list[str],
    gene_level_parallelism: int,
    force: bool,
    conditional_iterations: int = 4,
    joint_call_table_path: str = DEFAULT_JOINT_CALL_TABLE_PATH,
    frequency_table_path: str = DEFAULT_FREQUENCY_TABLE_PATH,
    vep_annotation_table_path: str = DEFAULT_VEP_ANNOTATION_TABLE_PATH,
    gencode_gtf_path: str = DEFAULT_GENCODE_GTF_PATH,
):
    """Run association script for all chromosomes and cell types"""

    if not any(
        input_files_prefix.startswith(prefix) for prefix in ('gs://', '/', 'https://')
    ):
        input_files_prefix = output_path(input_files_prefix)
    if not chromosomes:
        chromosomes = get_chromosomes(input_files_prefix)

    if not (isinstance(chromosomes, list) and len(chromosomes) > 0):
        raise ValueError('Must specify at least 1 chromosome as a list')

    if cell_types is None or len(cell_types) == 0:
        # not provided (ie: use all cell types)
        cell_types = get_cell_types_from(input_files_prefix)
        _logger.info(f'Found {len(cell_types)} cell types: {cell_types}')

    if not cell_types:
        raise ValueError(f'No cell types found at: {input_files_prefix}')

    # ideally this would come from metamist :(
    keys_tsv_path = os.path.join(input_files_prefix, 'OneK1K_CPG_IDs.tsv')

    outputs: dict[str, dict[str, Any]] = defaultdict(dict)

    # do the genotype_info stuff
    filtered_mt_path = output_path('genotype_table.mt', 'tmp')
    if AnyPath(filtered_mt_path).exists() and not force:
        _logger.info(f'Reusing filtered_mt path: {filtered_mt_path}')
    else:
        filter_mt_job = batch.new_python_job(f'filter_mt')
        copy_common_env(filter_mt_job)
        filtered_mt_path = filter_mt_job.call(
            filter_joint_call_mt,
            keys_path=keys_tsv_path,
            joint_mt_path=joint_call_table_path,
            frequency_table_path=frequency_table_path,
            vep_annotation_path=vep_annotation_table_path,
            output_location=filtered_mt_path,
            force=force,
        )

    for cell_type in cell_types:
        expression_tsv_path = os.path.join(
            input_files_prefix, 'expression_files', f'{cell_type}_expression.tsv'
        )
        covariates_tsv_path = os.path.join(
            input_files_prefix, 'covariates_files', f'{cell_type}_peer_factors_file.txt'
        )

        for chromosome in chromosomes:
            geneloc_tsv_path = os.path.join(
                input_files_prefix,
                'gene_location_files',
                f'GRCh38_geneloc_chr{chromosome}.tsv',
            )

            _genes = genes or get_genes_for_chromosome(
                expression_tsv_path=expression_tsv_path,
                geneloc_tsv_path=geneloc_tsv_path,
            )

            eqtl_outputs = generate_eqtl_spearman(
                batch=batch,
                # constants
                force=force,
                job_prefix=f'{cell_type}_chr{chromosome}_',
                cell_type=cell_type,
                chromosome=chromosome,
                output_prefix=output_path(f'{cell_type}/{chromosome}'),
                genes=_genes,
                gene_level_parallelism=gene_level_parallelism,
                # derived paths
                filtered_mt_path=filtered_mt_path,
                gencode_gtf_path=gencode_gtf_path,
                expression_tsv_path=expression_tsv_path,
                covariates_tsv_path=covariates_tsv_path,
                geneloc_tsv_path=geneloc_tsv_path,
            )

            conditional_outputs = generate_conditional_analysis(
                batch=batch,
                # constants
                force=force,
                job_prefix=f'{cell_type}_chr{chromosome}_conditional_',
                genes=_genes,
                gene_level_parallelism=gene_level_parallelism,
                cell_type=cell_type,
                chromosome=chromosome,
                filtered_matrix_table_path=filtered_mt_path,
                significant_snps_directory=eqtl_outputs['spearman_parquet_directory'],
                residuals_path=eqtl_outputs['residuals_csv_path'],
                iterations=conditional_iterations,
                output_prefix=output_path(f'{cell_type}/{chromosome}'),
            )

            outputs[cell_type][chromosome] = {
                **{'eqtl_' + k: v for k, v in eqtl_outputs.items()},
                **{'conditional_' + k: v for k, v in conditional_outputs.items()},
            }


# endregion DRIVER


if __name__ == '__main__':
    from_cli()  # pylint: disable=no-value-for-parameter
