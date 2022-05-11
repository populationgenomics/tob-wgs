#!/usr/bin/env python3

"""Perform conditional analysis on SNPs and expression residuals"""

import os
import hail as hl
import hailtop.batch as hb
import pandas as pd
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
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


DEFAULT_DRIVER_MEMORY = '4G'
DRIVER_IMAGE = os.getenv('CPG_DRIVER_IMAGE')
assert DRIVER_IMAGE

TOB_WGS = dataset_path('mt/v7.mt/')
FREQ_TABLE = dataset_path('joint-calling/v7/variant_qc/frequencies.ht/', 'analysis')


def get_number_of_scatters(residual_df, significant_snps_df):
    """get index of total number of genes"""

    # Identify the top eSNP for each eGene and assign remaining to df
    esnp1 = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .first()
        .reset_index()
    )
    gene_ids = esnp1['gene_symbol'][esnp1['gene_symbol'].isin(residual_df.columns)]

    return len(gene_ids)


def prepare_genotype_info(keys_path):

    init_batch()
    filtered_mt_path = output_path('genotype_table.ht', 'tmp')
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
        # filter out variants with MAF < 0.01
        ht = hl.read_table(FREQ_TABLE)
        mt = mt.annotate_rows(freq=ht[mt.row_key].freq)
        mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
        # add OneK1K IDs to genotype mt
        sampleid_keys = pd.read_csv(AnyPath(keys_path), sep='\t')
        print(f'Printing sampleid_keys: {sampleid_keys.head()}')
        genotype_samples = pd.DataFrame(mt.s.collect(), columns=['sampleid'])
        print(f'Printing genotype samples: {genotype_samples.head()}')
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
    """load genotype df and filter"""
    init_batch()
    mt = hl.read_matrix_table(filtered_mt_path)
    # only keep samples that are contained within the residuals df
    # this is important, since not all indivuduals have expression/residual
    # data (this varies by cell type)
    print(f'printing residual_df: {residual_df.head()}')
    print(f'printing residual_df columns: {residual_df.columns}')
    print(f'printing residual_df, dtype: {residual_df.dtypes}')
    samples_to_keep = set(residual_df.sampleid)
    print(f'printing samples_to_keep: {samples_to_keep}')
    set_to_keep = hl.literal(samples_to_keep)
    mt = mt.filter_cols(set_to_keep.contains(mt['onek1k_id']))
    print(f'printing mt: {mt.show()}')
    # Do this only on SNPs contained within gene_snp_df to save on
    # computational time
    snps_to_keep = set(gene_snp_test_df.snpid)
    print(f'printing snps_to_keep: {snps_to_keep}')
    sorted_snps = sorted(snps_to_keep)
    # only keep chromosome and position
    sorted_snp_positions = list(map(lambda x: x.split(':')[:2][1], sorted_snps))
    # convert all elements in list to int type
    sorted_snp_positions = [int(i) for i in sorted_snp_positions]
    print(f'printing sorted_snp_positions: {sorted_snp_positions}')
    # get first and last positions, with 1 added to last position (to make it inclusive)
    chromosome = gene_snp_test_df.snpid[0].split(':')[:1][0]
    first_and_last_snp = chromosome + ':' + str(sorted_snp_positions[0]) + '-' + str(sorted_snp_positions[-1]+1)
    print(f'printing first_and_last_snp: {first_and_last_snp}')
    # parse mt to region of interest
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(first_and_last_snp, reference_genome='GRCh38')])
    t = mt.entries()
    t = t.annotate(n_alt_alleles=t.GT.n_alt_alleles())
    t = t.key_by(contig=t.locus.contig, position=t.locus.position)
    t = t.select(t.alleles, t.onek1k_id, t.n_alt_alleles)
    t = t.annotate(
        snpid=hl.str(t.contig)
        + ':'
        + hl.str(t.position)
        + ':'
        + hl.str(t.alleles[0])
        + ':'
        + hl.str(t.alleles[1])
    )
    # Further reduce the table by only selecting SNPs needed
    set_to_keep = hl.literal(snps_to_keep)
    t = t.filter(set_to_keep.contains(t['snpid']))
    # only keep SNPs where all samples have an alt_allele value
    snps_to_remove = set(t.filter(hl.is_missing(t.n_alt_alleles)).snpid.collect())
    print(f'printing snps_to_remove: {snps_to_remove}')
    if len(snps_to_remove) > 0:
        t = t.filter(~hl.literal(snps_to_remove).contains(t.snpid))
        genotype_df = t.to_pandas(flatten=True)
        print(f'printing genotype_df: {genotype_df.head()}')
        genotype_df.rename({'onek1k_id': 'sampleid'}, axis=1, inplace=True)
        print(f'printing genotype_df: {genotype_df.head()}')
    else:
        genotype_df = t.to_pandas(flatten=True)
        print(f'printing genotype_df: {genotype_df.head()}')
        genotype_df.rename({'onek1k_id': 'sampleid'}, axis=1, inplace=True)
        print(f'printing genotype_df: {genotype_df.head()}')

    print(f'printing genotype_df: {genotype_df.head()}')
    return genotype_df


def calculate_residual_df(residual_df, significant_snps_df, filtered_mt_path):
    """calculate residuals for gene list"""
    
    print(f'printing residual_df within loop: {residual_df.head()}')
    print(f'printing residual_df within loop, dtype: {residual_df.dtypes}')

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
    print(f'printing residual_df within loop, dtype: {residual_df.dtypes}')
    print(f'printing residual_df within loop: {residual_df.head()}')

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
        snp = esnp1.snpid[esnp1.gene_symbol == gene].to_string(index=False)
        snp_genotype = genotype_df[genotype_df.snpid == snp][['sampleid', 'n_alt_alleles']]

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
    print(f'printing adjusted_residual_mat within loop: {adjusted_residual_mat.head()}')
    print(f'printing adjusted_residual_mat within loop, dtype: {adjusted_residual_mat.dtypes}')
    adjusted_residual_mat.columns = gene_ids
    adjusted_residual_mat.insert(loc=0, column='sampleid', value=genotype_df.sampleid)
    print(f'printing adjusted_residual_mat within loop: {adjusted_residual_mat.head()}')

    # call adjusted_residual_mat.sampleid - should fail
    return adjusted_residual_mat


# Run Spearman rank in parallel by sending genes in batches
def run_computation_in_scatter(
    iteration,  # pylint: disable=redefined-outer-name
    idx,
    residual_df,
    significant_snps_df,
    filtered_mt_path
):
    """Run genes in scatter"""

    print(f'iteration = {iteration+2}')
    print(f'idx = {idx}')

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
    gene_snp_test_df = esnps_to_test[['snpid', 'gene_symbol', 'gene_id']]
    gene_snp_test_df = gene_snp_test_df[
        gene_snp_test_df['gene_symbol'] == gene_ids.iloc[idx]
    ]
    print(f'printing gene_snp_test_df: {gene_snp_test_df.head()}')
    # Subset genotype file for the significant SNPs
    genotype_df = get_genotype_df(
        filtered_mt_path, residual_df, gene_snp_test_df
    )
    print(f'Printing genotype_df {genotype_df.head()}')

    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_symbol = df.gene_symbol
        gene_id = df.gene_id
        snp = df.snpid
        print(snp)
        gt = genotype_df[genotype_df.snpid == snp][['sampleid', 'n_alt_alleles']]
        print(f'Printing gt {gt.head()}')

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
        'snpid',
        'spearmans_rho',
        'p_value',
    ]
    # remove any NA values. Remove this line when run on main dataset
    adjusted_spearman_df = adjusted_spearman_df.dropna(axis=0, how='any')
    # add in global position and round
    locus = adjusted_spearman_df.snpid.str.split(':', expand=True)[[0, 1]].agg(
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
    adjusted_spearman_df['round'] = iteration + 2
    print(f'adjusted_spearman_df finished')
    print(f'adjusted_spearman_df.head(): {adjusted_spearman_df.head()}')
    # init_batch()
    t = hl.Table.from_pandas(adjusted_spearman_df)
    print(f't.show(): {t.show()}')
    t = t.annotate(global_bp=hl.locus(t.chrom, hl.int32(t.bp)).global_position())
    t = t.annotate(locus=hl.locus(t.chrom, hl.int32(t.bp)))
    # get alleles
    mt = hl.read_matrix_table(filtered_mt_path).key_rows_by('locus')
    t = t.key_by('locus')
    t = t.annotate(
        a1=mt.rows()[t.locus].alleles[0],
        a2=mt.rows()[t.locus].alleles[1],
    )
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
    adjusted_spearman_df = t.to_pandas()

    # set variables for next iteration of loop
    significant_snps_df = adjusted_spearman_df

    return significant_snps_df


def merge_significant_snps_dfs(*df_list):
    """
    Merge list of list of sig_snps dataframes
    """

    merged_sig_snps = pd.concat(df_list)
    pvalues = merged_sig_snps['p_value']
    fdr_values = pd.DataFrame(list(multi.fdrcorrection(pvalues))).iloc[1]
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
    '--significant-snps', required=True, help='A space separated list of SNPs'
)
@click.option('--residuals', required=True, help='A CSV of gene residuals')
@click.option(
    '--keys',
    required=True,
    help='A TSV of sample ids to convert external to internal IDs',
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
    Creates a Hail Batch pipeline for calculating EQTL using {iterations} iterations,
    scattered across the number of genes. Note, {iterations} iterations are run, however
    iterations start at 2, since round 1 is completed in `generate_eqtl_spearan.py`.
    """
    dataset = os.getenv('CPG_DATASET')
    backend = hb.ServiceBackend(billing_project=dataset, remote_tmpdir=remote_tmpdir())
    batch = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

    residual_df = pd.read_csv(AnyPath(residuals))
    significant_snps_df = pd.read_csv(AnyPath(
        significant_snps), sep=' ', skipinitialspace=True
    )
    
    if test_subset_genes:
        n_genes = test_subset_genes
    else:
        print('Loaded data to prepare workflow')
        # test with 5 genes
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
    for iteration in range(iterations):

        calc_resid_df_job = batch.new_python_job(
            f'calculate-resid-df-iter-{iteration+2}'
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
            os.path.join(output_prefix, f'round{iteration+2}_residual_results.csv'),
        )

        sig_snps_dfs = []
        for gene_idx in range(n_genes):
            j = batch.new_python_job(name=f'process_iter_{iteration+2}_job_{gene_idx}')
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
                filtered_mt_path
            )
            sig_snps_dfs.append(gene_result)

        merge_job = batch.new_python_job(name='merge_scatters')
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
            output_prefix, f'esnp_round{iteration+2}_table.csv'
        )
        batch.write_output(sig_snps_as_string.as_str(), sig_snps_output_path)

    batch.run(wait=False)


if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    main()
