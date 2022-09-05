#!/usr/bin/env python3

"""Run Spearman rank correlation on SNPs and expression residuals"""

import os
import click
import hail as hl
import hailtop.batch as hb
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from cpg_utils.hail_batch import (
    output_path,
    copy_common_env,
    init_batch,
    remote_tmpdir,
    dataset_path,
)
from cpg_utils.config import get_config
from cloudpathlib import AnyPath

from constants import (
    DEFAULT_JOINT_CALL_TABLE_PATH,
    DEFAULT_FREQUENCY_TABLE_PATH,
    DEFAULT_VEP_ANNOTATION_TABLE_PATH,
    DEFAULT_GENCODE_GTF_PATH,
    MULTIPY_IMAGE,
)

DEFAULT_DRIVER_MEMORY = '4G'


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


def get_number_of_scatters(*, expression_tsv_path, geneloc_tsv_path):
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
    # TODO: remove this
    return 296

    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep='\t')
    geneloc_df = pd.read_csv(AnyPath(geneloc_tsv_path), sep='\t')

    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = list(expression_df.columns.values)[1:]
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(gene_ids)]

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
    log_expression_df.insert(loc=0, column='sampleid', value=sample_ids)
    # remove sc outlier samples
    log_expression_df = remove_sc_outliers(log_expression_df)

    return log_expression_df


def calculate_log_cpm(expression_df):
    expression_df = filter_lowly_expressed_genes(expression_df)
    # remove sampleid column and get log expression
    # this can only be done on integers
    expression_df = expression_df.iloc[:, 1:]
    cpm_df = expression_df.apply(lambda x: (x / sum(x)) * 1000000, axis=0)
    log_cpm = np.log(cpm_df + 1)

    return log_cpm


def generate_log_cpm_output(
    *,
    output_prefix: str,
    cell_type: str,
    expression_tsv_path: str,
    gencode_gtf_path: str,
):
    """Calculate log cpm for each cell type and chromosome

    Input:
    expression_df: a data frame with samples as rows and genes as columns,
    containing normalised expression values (i.e., the average number of molecules
    for each gene detected in each person).

    Returns:
    Counts per million mapped reads
    """

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
    # convert int to str in order to avoid "int() argument must be a string, a bytes-like object or a number, not 'NoneType'"
    gtf = gtf.annotate(frame=hl.str(gtf.frame))
    gtf = gtf.to_pandas()
    data_summary['ensembl_ids'] = data_summary.merge(
        gtf.drop_duplicates('gene_name'), how='left', on='gene_name'
    ).gene_id

    # Save file
    data_summary_path = AnyPath(output_prefix) / 'gene_expression.parquet'
    with data_summary_path.open('wb') as fp:
        data_summary.to_parquet(fp)

    return data_summary_path.as_uri()


def calculate_residuals(expression_tsv_path, covariates_tsv_path, output_prefix):
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
    residual_path = AnyPath(output_prefix) / 'log_residuals.csv'
    with residual_path.open('w') as fp:
        residual_df.to_csv(fp, index=False)

    return residual_path.as_uri()


# Run Spearman rank in parallel by sending genes in batches
def run_spearman_correlation_scatter(
    idx,
    *,
    expression_tsv_path,
    geneloc_tsv_path,
    residuals_csv_path,
    filtered_mt_path,
    cell_type,
    parquet_output_dir: str,
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

    # import multipy here to avoid issues with driver image updates
    from multipy.fdr import qvalue

    # calculate log expression values
    residuals_df = pd.read_csv(AnyPath(residuals_csv_path))
    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep='\t')
    log_expression_df = get_log_expression(expression_df)

    # TODO: why calculate the sliding window around each gene,
    #   then filter for the specific gene, why not get the gene first?
    # Get 1Mb sliding window around each gene
    geneloc_df = pd.read_csv(AnyPath(geneloc_tsv_path), sep='\t')
    gene_ids = list(log_expression_df.columns.values)[1:]
    geneloc_df = geneloc_df[geneloc_df.gene_name.isin(gene_ids)]
    geneloc_df = geneloc_df.assign(left=geneloc_df.start - 1000000)
    geneloc_df = geneloc_df.assign(right=geneloc_df.end + 1000000)

    # perform correlation in chunks by gene
    gene_info = geneloc_df.iloc[idx]
    gene = gene_info.gene_name

    # checkpoint
    spearman_parquet_path = (
        AnyPath(parquet_output_dir) / f'correlation_results_{gene}.parquet'
    )
    if spearman_parquet_path.exists() and not force:
        return spearman_parquet_path

    chromosome = gene_info.chr
    # get all SNPs which are within 1Mb of each gene
    init_batch(driver_cores=8, worker_cores=8)
    mt = hl.read_matrix_table(filtered_mt_path)
    # only keep samples that are contained within the residuals df
    # this is important, since not all individuals have expression/residual
    # data (this varies by cell type)
    # this also filters out any sc sample outliers
    samples_to_keep = hl.literal(list(residuals_df.sampleid))
    mt = mt.filter_cols(samples_to_keep.contains(mt['onek1k_id']))

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

    association_effect_data = get_association_effect_data(gene)

    # Save file
    # TODO: hmm, this feels a little hacky, is it actually tmp, can we write
    #   to a tmp bucket?
    # tmp_dir = output_prefix.replace(
    #     output_prefix.split('/')[2], output_prefix.split('/')[2] + '-tmp'
    # )
    # edit: short term I've replaced it
    path = AnyPath(output_path(f'eqtl_effect_{gene}.parquet', 'tmp'))
    with path.open('wb') as fp:
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
        return (gene_symbol, gene_id, alleles, snp, coef, p)

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
    with spearman_parquet_path.open('wb') as fp:
        spearman_df.to_parquet(fp)

    return spearman_parquet_path


# Create click command line to enter dependency files
@click.command()
@click.option(
    '--expression',
    required=True,
    help='A TSV of normalised expression values, \
        where values represent the average number of molecules for each gene. Gene information \
        is contained in columns and sample information in rows.',
)
@click.option(
    '--geneloc',
    required=True,
    help='A TSV of start and end positions for each gene. Each row in \
        this df contains a gene, with ENSEMBL gene id, gene symbol, chromosome, start and end \
        positions, and strand information as columns.',
)
@click.option(
    '--covariates',
    required=True,
    help='A TSV of covariates to calculate residuals. Sample ids \
        are contained in rows and covariate information in columns. As a minimum, covariate \
        information should include PCA scores for the first 4 PCs, the top 2 PEER factors, \
        sex, and age.',
)
@click.option(
    '--filtered-mt',
    required=True,
    help='Filtered MT, default location is: output_path("genotype_table.mt", "tmp")',
)
def from_cli(expression, geneloc, covariates, filtered_mt: str):

    cell_type = expression.split('/')[-1].split('_expression')[0]
    chromosome = geneloc.split('_')[-1].removesuffix('.tsv')

    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name=f'Generate EQTL spearman ({cell_type}:',
        backend=backend,
        default_python_image=MULTIPY_IMAGE,
    )

    main(
        batch=batch,
        job_prefix='',
        output_prefix=output_path(f'{cell_type}/{chromosome}'),
        chromosome=chromosome,
        cell_type=cell_type,
        # input paths
        filtered_mt_path=filtered_mt,
        expression_tsv_path=expression,
        geneloc_tsv_path=geneloc,
        covariates_tsv_path=covariates,
        gencode_gtf_path=DEFAULT_GENCODE_GTF_PATH,
    )

    batch.run(wait=False)


def main(
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
    cell_type: str,
    chromosome: str,
    force: bool = False,
):
    """
    Creates a Hail Batch pipeline for calculating EQTLs
    """

    n_genes_in_scatter = get_number_of_scatters(
        expression_tsv_path=expression_tsv_path, geneloc_tsv_path=geneloc_tsv_path
    )
    if n_genes_in_scatter == 0:
        # only run batch jobs if there are test genes in the chromosome
        raise ValueError('No genes in get_number_of_scatters()')

    calculate_residuals_job = batch.new_python_job(f'{job_prefix}calculate-residuals')
    copy_common_env(calculate_residuals_job)
    residuals_csv_path = calculate_residuals_job.call(
        calculate_residuals,
        output_prefix=output_prefix,
        expression_tsv_path=expression_tsv_path,
        covariates_tsv_path=covariates_tsv_path,
    )
    # calculate and save log cpm values
    generate_log_cpm_output_job = batch.new_python_job(
        f'{job_prefix}generate_log_cpm_output'
    )
    generate_log_cpm_output_job.memory('8Gi')
    copy_common_env(generate_log_cpm_output_job)
    generate_log_cpm_output_job.call(
        generate_log_cpm_output,
        cell_type=cell_type,
        output_prefix=output_prefix,
        expression_tsv_path=expression_tsv_path,
        gencode_gtf_path=gencode_gtf_path,
    )

    spearman_output_directory = os.path.join(output_prefix, 'parquet')
    spearman_output_dependencies = []

    for idx in range(n_genes_in_scatter):
        j = batch.new_python_job(
            name=f'{job_prefix}calculate_spearman_correlation_{idx}'
        )
        j.depends_on(calculate_residuals_job)
        j.cpu(2)
        j.memory('8Gi')
        j.storage('2Gi')
        j.image(MULTIPY_IMAGE)
        copy_common_env(j)
        j.call(
            run_spearman_correlation_scatter,
            idx=idx,
            expression_tsv_path=expression_tsv_path,
            geneloc_tsv_path=geneloc_tsv_path,
            residuals_csv_path=residuals_csv_path,
            filtered_mt_path=filtered_mt_path,
            celltype=cell_type,
            parquet_output_dir=spearman_output_directory,
            force=force,
        )
        # we'll track job dependencies with an explicit sink
        spearman_output_dependencies.append(j)

    # create a fake job to avoid having to carry dependencies outside this script
    def return_value(value):
        return value

    fake_sink = batch.new_python_job(f'{job_prefix}sink')
    fake_sink.depends_on(*spearman_output_dependencies)
    sinked_output_dir = fake_sink.call(return_value, spearman_output_directory)

    return {
        'spearman_parquet_directory': sinked_output_dir,
        'residuals_csv_path': residuals_csv_path,
    }


if __name__ == '__main__':
    from_cli()  # pylint: disable=no-value-for-parameter
