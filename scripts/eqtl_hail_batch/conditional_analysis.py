#!/usr/bin/env python3

"""Perform conditional analysis on SNPs and expression residuals"""

import os
import hail as hl
import hailtop.batch as hb
import pandas as pd
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
from cpg_utils.config import get_config
from cloudpathlib import AnyPath
import click


DEFAULT_DRIVER_MEMORY = '4G'
MULTIPY_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16'

TOB_WGS = dataset_path('mt/v7.mt/')
FREQ_TABLE = dataset_path('joint-calling/v7/variant_qc/frequencies.ht/', 'analysis')


def get_number_of_scatters(residual_df: pd.DataFrame, significant_snps_df: pd.DataFrame) -> int:
    """get index of total number of genes
    
    Input:
    significant_snps_df: a dataframe with significance results (p-value and FDR) from the Spearman rank 
    correlation calculated in `conditional_analysis.py`. Each row contains information for a SNP, 
    and columns contain information on significance levels, spearmans rho, and gene information.

    Returns:
    The number of genes (as an int) returned after sorting each row of the significant_snps_df
    and taking the top SNP for each gene. This integer number gets fed into the number of 
    scatters to run.
    """

    # Identify the top eSNP for each eGene
    esnp1 = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .first()
        .reset_index()
    )
    gene_ids = esnp1['gene_symbol'][esnp1['gene_symbol'].isin(residual_df.columns)]

    return len(gene_ids)


def prepare_genotype_info(keys_path):
    """Filter hail matrix table

    Input:
    keys_path: path to a tsv file with information on
    OneK1K amd CPG IDs (columns) for each sample (rows).

    Returns:
    Path to a hail matrix table, with rows (alleles) filtered on the following requirements:
    1) biallelic, 2) meets VQSR filters, 3) gene quality score higher than 20,
    4) call rate of 0.8, and 5) variants with MAF <= 0.01. Columns (samples) are filtered
    on the basis of having rna-seq expression data, i.e., within the filtered log_expression_df
    """

    init_batch()
    filtered_mt_path = output_path('genotype_table.mt', 'tmp')
    if not hl.hadoop_exists(filtered_mt_path):
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
        mt.write(filtered_mt_path)

    return filtered_mt_path


def get_genotype_df(filtered_mt_path, residual_df, gene_snp_test_df):
    """load genotype df and filter
    Input:
    filtered_mt_path: the path to a filtered matrix table generated using prepare_genotype_info(), 
    (and filtering requirements outlined in the same function).
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
    mt = hl.read_matrix_table(filtered_mt_path)
    # only keep samples that are contained within the residuals df
    # this is important, since not all indivuduals have expression/residual
    # data (this varies by cell type)
    samples_to_keep = set(residual_df.sampleid)
    set_to_keep = hl.literal(samples_to_keep)
    mt = mt.filter_cols(set_to_keep.contains(mt['onek1k_id']))
    # Do this only on SNPs contained within gene_snp_df to save on
    # computational time
    snps_to_keep = set('chr' + gene_snp_test_df.snp_id)
    sorted_snps = sorted(snps_to_keep)
    sorted_snp_positions = list(map(lambda x: x.split(':')[:2][1], sorted_snps))
    sorted_snp_positions = [int(i) for i in sorted_snp_positions]
    # get first and last positions, with 1 added to last position (to make it inclusive)
    chromosome = 'chr' + gene_snp_test_df.snp_id[0].split(':')[:1][0]
    first_and_last_snp = chromosome + ':' + str(sorted_snp_positions[0]) + '-' + str(sorted_snp_positions[-1]+1)
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(first_and_last_snp, reference_genome='GRCh38')])
    t = mt.entries()
    t = t.annotate(n_alt_alleles=t.GT.n_alt_alleles())
    t = t.key_by(contig=t.locus.contig, position=t.locus.position)
    t = t.select(t.alleles, t.n_alt_alleles, sampleid=t.onek1k_id)
    t = t.annotate(snp_id=hl.str(':').join(list(map(hl.str, [t.contig, t.position, t.alleles[0], t.alleles[1]]))))
    # Further reduce the table by only selecting SNPs needed
    set_to_keep = hl.literal(snps_to_keep)
    t = t.filter(set_to_keep.contains(t['snp_id']))
    # only keep SNPs where all samples have an alt_allele value
    snps_to_remove = set(t.filter(hl.is_missing(t.n_alt_alleles)).snp_id.collect())
    if len(snps_to_remove) > 0:
        t = t.filter(~hl.literal(snps_to_remove).contains(t.snp_id))

    genotype_df = t.to_pandas(flatten=True)
    genotype_df.rename({'onek1k_id': 'sampleid'}, axis=1, inplace=True)
    # chromosome must be turned into a number solely
    # note, 'chr' + chromosome number was necessary in the previous step, 
    # as 'chr' indicates GrCh38 format
    genotype_df['contig'] = genotype_df.contig.str.split('chr', expand=True)[1]
    genotype_df['snp_id'] = genotype_df.snp_id.str.split('chr', expand=True)[1]

    return genotype_df


def calculate_residual_df(residual_df, significant_snps_df, filtered_mt_path):
    """calculate residuals for gene list

    Input:
    residual_df: a dataframe of expression residuals, with genes as columns and samples
    as rows. 
    significant_snps_df: a dataframe with significance results (p-value and FDR) from the Spearman rank 
    correlation calculated in `conditional_analysis.py`. Each row contains information for a SNP, 
    and columns contain information on significance levels, spearmans rho, and gene information.
    filtered_mt_path: the path to a filtered matrix table generated using prepare_genotype_info(), 
    (and filtering requirements outlined in the same function).

    Returns: a dataframe of expression residuals, with genes as columns and samples
    as rows, which have been conditioned on the lead eSNP. 
    """
    
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

    # Subset residuals for the genes to be tested
    gene_ids = esnp1['gene_symbol'][esnp1['gene_symbol'].isin(residual_df.columns)].to_list()
    # save sampleids before filtering redidual_df
    sample_ids = residual_df.sampleid
    residual_df = residual_df[gene_ids]
    # reassign sample ids
    residual_df = residual_df.assign(sampleid=sample_ids.to_list())

    # Subset genotype file for the significant SNPs
    genotype_df = get_genotype_df(
        filtered_mt_path, residual_df, esnp1
    )

    # Find residuals after adjustment of lead SNP
    def calculate_adjusted_residuals(gene_id):
        gene = gene_id
        # select gene to regress
        exprs_val = residual_df[['sampleid', gene]]
        # select SNP to add
        snp = esnp1.snp_id[esnp1.gene_symbol == gene].to_string(index=False)
        snp_genotype = genotype_df[genotype_df.snp_id == snp][['sampleid', 'n_alt_alleles']]

        # Create a test df by adding covariates
        test_df = exprs_val.merge(snp_genotype, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'expression', 'genotype']

        y, x = dmatrices('expression ~ genotype', test_df)
        model = sm.OLS(y, x)
        residuals = list(model.fit().resid)
        return residuals

    adjusted_residual_mat = pd.DataFrame(
        list(map(calculate_adjusted_residuals, gene_ids))
    ).T
    adjusted_residual_mat.columns = gene_ids
    adjusted_residual_mat.insert(loc=0, column='sampleid', value=genotype_df.sampleid)

    return adjusted_residual_mat


# Run Spearman rank in parallel by sending genes in batches
def run_computation_in_scatter(
    iteration,  # pylint: disable=redefined-outer-name
    idx,
    residual_df,
    significant_snps_df,
    filtered_mt_path,
    celltype,
    output_prefix
):
    """Run genes in scatter
    
    Input:
    iteration: synonymous with round. 4 iterations are performed by default for the conditional
    analysis, which results in 5 total rounds with the addition of the first analysis. 
    idx: the index of each gene with a top eSNP. Note: a 1Mb region is taken 
    up and downstream of each gene.
    residual_df: residual dataframe, calculated in calculate_residual_df(). This is run for
    each round/iteration of conditional analysis. 
    filtered_mt_path: the path to a filtered matrix table generated using prepare_genotype_info(), 
    (and filtering requirements outlined in the same function).

    Returns:
    A dataframe with significance results (p-value and FDR) from the Spearman rank correlation. Each
    row contains information for a SNP, and columns contain information on significance levels, 
    spearmans rho, and gene information.
    """

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
    esnp1_path = AnyPath(output_prefix) / f'conditioned_esnps_{iteration}.tsv'
    with esnp1_path.open('w') as fp:
        esnp1.to_csv(fp, index=False)

    # Remaning eSNPs to test
    esnps_to_test = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .apply(lambda group: group.iloc[1:, 1:])
        .reset_index()
    )

    # for each gene, get esnps_to_test
    gene_ids = esnp1['gene_symbol'][esnp1['gene_symbol'].isin(residual_df.columns)]
    esnps_to_test = esnps_to_test[esnps_to_test.gene_symbol.isin(residual_df.columns)]
    gene_snp_test_df = esnps_to_test[['snp_id', 'gene_symbol', 'gene_id']]
    gene_snp_test_df = gene_snp_test_df[
        gene_snp_test_df['gene_symbol'] == gene_ids.iloc[idx]
    ]
    # Subset genotype file for the significant SNPs
    genotype_df = get_genotype_df(
        filtered_mt_path, residual_df, gene_snp_test_df
    )

    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_symbol = df.gene_symbol
        gene_id = df.gene_id
        snp = df.snp_id
        gt = genotype_df[genotype_df.snp_id == snp][['sampleid', 'n_alt_alleles']]

        res_val = residual_df[['sampleid', gene_symbol]]
        test_df = res_val.merge(gt, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        spearmans_rho, p = spearmanr(
            test_df['SNP'], test_df['residual'], nan_policy='omit'
        )
        return (gene_symbol, gene_id, snp, spearmans_rho, p)

    # calculate spearman correlation
    adjusted_spearman_df = pd.DataFrame(
        list(gene_snp_test_df.apply(spearman_correlation, axis=1))
    )
    adjusted_spearman_df.columns = [
        'gene_symbol',
        'gene_id',
        'snp_id',
        'spearmans_rho',
        'p_value',
    ]
    # remove any NA values. Remove this line when run on main dataset
    adjusted_spearman_df = adjusted_spearman_df.dropna(axis=0, how='any')
    # add in global position and round
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
    t = hl.Table.from_pandas(adjusted_spearman_df)
    t = t.annotate(global_bp=hl.locus(t.chrom, hl.int32(t.bp)).global_position())
    t = t.annotate(locus=hl.locus(t.chrom, hl.int32(t.bp)))
    # get alleles
    mt = hl.read_matrix_table(filtered_mt_path).key_rows_by('locus')
    t = t.key_by('locus')
    t = t.annotate(
        a1=mt.rows()[t.locus].alleles[0],
        a2=mt.rows()[t.locus].alleles[1],
    )
    # turn back into pandas df and add additional information
    # for front-end analysis
    adjusted_spearman_df = t.to_pandas()
    # add celltype id
    celltype_id = celltype.lower()
    adjusted_spearman_df['cell_type_id'] = celltype_id
    adjusted_spearman_df['is_esnp'] = adjusted_spearman_df.snp_id.isin(esnp1.snp_id)
    # chromosome must be turned into a number solely
    # note, 'chr' + chromosome number was necessary in the previous step, 
    # as 'chr' indicates GrCh38 format
    adjusted_spearman_df['chrom'] = adjusted_spearman_df.chrom.str.split('chr', expand=True)[1]
    # add association ID annotation after adding in alleles, a1, and a2
    adjusted_spearman_df['association_id'] = adjusted_spearman_df.apply(lambda x: ':'.join(x[['chrom', 'bp', 'a1', 'a2', 'gene_symbol', 'cell_type_id', 'round']]), axis=1)
    adjusted_spearman_df['variant_id'] = adjusted_spearman_df.apply(lambda x: ':'.join(x[['chrom', 'bp', 'a2']]), axis=1)
    adjusted_spearman_df['snp_id'] = adjusted_spearman_df.apply(lambda x: ':'.join(x[['chrom', 'bp', 'a1', 'a2']]), axis=1)

    # set variables for next iteration of loop
    significant_snps_df = adjusted_spearman_df

    return significant_snps_df


def merge_significant_snps_dfs(*df_list):
    """
    Merge list of list of sig_snps dataframes
    """

    # This PythonJob is run in the multipy container, 
    # do the import here so it's not run in the driver container
    from multipy.fdr import qvalue 

    merged_sig_snps = pd.concat(df_list)
    pvalues = merged_sig_snps['p_value']
    _, qvals = qvalue(pvalues)
    fdr_values = pd.DataFrame(list(qvals)).iloc[1]
    merged_sig_snps = merged_sig_snps.assign(fdr=fdr_values)
    merged_sig_snps['fdr'] = merged_sig_snps.fdr.astype(float)
    merged_sig_snps.append(merged_sig_snps)

    return merged_sig_snps


def convert_dataframe_to_text(dataframe):
    """
    convert to string for writing
    """
    return dataframe.to_string()


# Create click command line to enter dependency files
@click.command()
@click.option(
    '--significant-snps', required=True, help='A space separated list of SNPs (as rows), \
        and significance levels, spearmans rho, and gene information as columns.'
)
@click.option('--residuals', required=True, help='A CSV of gene residuals, with genes \
    as columns and samples as rows.')
@click.option(
    '--keys',
    required=True,
    help='A TSV of sample ids to convert external to internal IDs. Rows contain \
        sample ids, while columns contain OneK1K IDs and CPG IDs.',
) 
@click.option(
    '--output-prefix',
    required=True,
    help='A path prefix of where to output files, eg: gs://MyBucket/output-folder/',
)
@click.option(
    '--iterations', type=int, default=4, help='Number of iterations to perform'
)
@click.option(
    '--test-subset-genes',  # pylint: disable=too-many-locals
    type=int,
    help='Test with {test_subset_genes} genes, often = 5.',
)
def main(
    significant_snps: str,
    residuals,
    keys,
    output_prefix: str,
    iterations=4,
    test_subset_genes=None,
):
    """
    Creates a Hail Batch pipeline for calculating eQTLs using {iterations} iterations,
    scattered across the number of genes. Note, iteration 1 is completed in `generate_eqtl_spearan.py`.
    """
    backend = hb.ServiceBackend(billing_project=get_config()['hail']['billing_project'], remote_tmpdir=remote_tmpdir())
    batch = hb.Batch(name='eQTL', backend=backend, default_python_image=get_config()['workflow']['driver_image'])

    # get cell type to feed into run_computation_in_scatter
    celltype = output_prefix.split('scrna-seq/')[-1].split('/')[0]

    residual_df = pd.read_csv(AnyPath(residuals))
    significant_snps_df = pd.read_csv(AnyPath(
        significant_snps), sep=' ', skipinitialspace=True
    )
    
    if test_subset_genes:
        n_genes = test_subset_genes
    else:
        n_genes = test_subset_genes or get_number_of_scatters(
            residual_df, significant_snps_df
        )

    filter_mt_job = batch.new_python_job('filter_mt')
    copy_common_env(filter_mt_job)
    filtered_mt_path = filter_mt_job.call(
        prepare_genotype_info, keys_path=keys
    )

    previous_sig_snps_result = significant_snps_df  # pylint: disable=invalid-name
    previous_residual_result = residual_df  # pylint: disable=invalid-name
    # Perform conditional analysis for n iterations (specified in click interface)
    # note, iteration 1 is performed in generate_eqtl_spearman.py, which requires starting
    # the iteration at 2 (from a 0 index)
    for iteration in range(2, iterations + 2):

        calc_resid_df_job = batch.new_python_job(
            f'calculate-resid-df-iter-{iteration}'
        )
        calc_resid_df_job.cpu(2)
        calc_resid_df_job.memory('8Gi')
        calc_resid_df_job.storage('2Gi')
        copy_common_env(calc_resid_df_job)
        previous_residual_result = calc_resid_df_job.call(
            calculate_residual_df,
            previous_residual_result,
            previous_sig_snps_result,
            filtered_mt_path
        )

        # convert residual df to string for output
        residual_as_str = calc_resid_df_job.call(
            convert_dataframe_to_text, previous_residual_result
        )
        # output residual df for each iteration
        batch.write_output(
            residual_as_str.as_str(),
            os.path.join(output_prefix, f'round{iteration}_residual_results.csv'),
        )

        sig_snps_dfs = []
        for gene_idx in range(n_genes):
            j = batch.new_python_job(name=f'process_iter_{iteration}_job_{gene_idx}')
            j.cpu(2)
            j.memory('8Gi')
            j.storage('2Gi')
            copy_common_env(j)
            gene_result: hb.resource.PythonResult = j.call(
                run_computation_in_scatter,
                iteration,
                gene_idx,
                previous_residual_result,
                previous_sig_snps_result,
                filtered_mt_path,
                celltype,
                output_prefix
            )
            sig_snps_dfs.append(gene_result)

        merge_job = batch.new_python_job(name='merge_scatters')
        merge_job.image(MULTIPY_IMAGE)
        merge_job.cpu(2)
        merge_job.memory('8Gi')
        merge_job.storage('2Gi')
        previous_sig_snps_result = merge_job.call(
            merge_significant_snps_dfs, *sig_snps_dfs
        )

        # convert sig snps to string for output
        sig_snps_as_string = merge_job.call(
            convert_dataframe_to_text, previous_sig_snps_result
        )
        # output sig snps for each iteration
        sig_snps_output_path = os.path.join(
            output_prefix, f'esnp_round{iteration}_table.csv'
        )
        batch.write_output(sig_snps_as_string.as_str(), sig_snps_output_path)

    batch.run(wait=False)


if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    main()
