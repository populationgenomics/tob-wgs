#!/usr/bin/env python3

"""Perform conditional analysis on SNPs and expression residuals"""
from typing import Any

import os
import hail as hl
import hailtop.batch as hb
import pandas as pd
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
FILTERED_MT = output_path('genotype_table.mt', 'tmp')


def get_number_of_scatters(*, residuals_path: str, significant_snps_path: str) -> int:
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
    if any(
        isinstance(p, hb.PythonResult) for p in (residuals_path, significant_snps_path)
    ):
        raise ValueError(
            'Residuals dir and significant snps MUST be a string, but '
            'received a hb.PythonResult, this means it\'s an intermediate '
            'result and cannot be operated on.'
        )

    residual_df = pd.read_csv(AnyPath(residuals_path))
    significant_snps_df = pd.read_parquet(AnyPath(significant_snps_path))

    # Identify the top eSNP for each eGene
    esnp1 = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .first()
        .reset_index()
    )
    gene_ids = esnp1['gene_symbol'][esnp1['gene_symbol'].isin(residual_df.columns)]

    return len(gene_ids)


def get_genotype_df(residual_df, gene_snp_test_df):
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
    mt = hl.read_matrix_table(FILTERED_MT)
    # only keep samples that are contained within the residuals df
    # this is important, since not all indivuduals have expression/residual
    # data (this varies by cell type)
    samples_to_keep = hl.set(list(residual_df.sampleid))
    mt = mt.filter_cols(samples_to_keep.contains(mt['onek1k_id']))

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


def calculate_residual_df(
    residual_path: str,
    significant_snps_path: str,
    output_path: str,
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

    if AnyPath(output_path).exists() and not force:
        return output_path

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
    genotype_df = get_genotype_df(residual_df, esnp1)

    # Find residuals after adjustment of lead SNP
    def calculate_adjusted_residuals(gene_id):
        gene = gene_id
        # select gene to regress
        exprs_val = residual_df[['sampleid', gene]]
        # select SNP to add
        snp = esnp1.snp_id[esnp1.gene_symbol == gene].to_string(index=False)
        snp_genotype = genotype_df[genotype_df.snp_id == snp][
            ['sampleid', 'n_alt_alleles']
        ]

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

    adjusted_residual_mat.to_csv(output_path, index=False)
    return output_path


# Run Spearman rank in parallel by sending genes in batches
def run_computation_in_scatter(
    iteration: int,
    idx: int,
    residual_path: str,
    significant_snps_path: str,
    cell_type: str,
    round_outputdir: str,
    sigsnps_outputdir: str,
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

    output_path = os.path.join(sigsnps_outputdir, f'sig-snps-{idx}.parquet')

    if AnyPath(output_path).exists() and not force:
        return output_path

    # import these here to avoid polluting the global namespace
    from multipy.fdr import qvalue
    from scipy.stats import spearmanr

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
    # save esnp1 for front-end use on which SNPs have been conditioned on
    esnp1_path = AnyPath(round_outputdir) / f'conditioned_esnps_{iteration}.tsv'
    with esnp1_path.open('w') as fp:
        esnp1.to_csv(fp, index=False)

    # Remaning eSNPs to test
    esnps_to_test = (
        significant_snps_df.sort_values(['gene_symbol', 'fdr'], ascending=True)
        .groupby('gene_symbol')
        .apply(lambda group: group.iloc[1:, 1:])
        .reset_index()
    )
    # add in SNP_ID column
    esnps_to_test['snp_id'] = esnps_to_test[['chrom', 'bp', 'a1', 'a2']].apply(
        lambda row: ':'.join(row.values.astype(str)), axis=1
    )

    # for each gene, get esnps_to_test
    gene_ids = esnp1['gene_symbol'][esnp1['gene_symbol'].isin(residual_df.columns)]
    esnps_to_test = esnps_to_test[esnps_to_test.gene_symbol.isin(residual_df.columns)]
    gene_snp_test_df = esnps_to_test[['snp_id', 'gene_symbol', 'gene_id', 'a1', 'a2']]
    gene_snp_test_df = gene_snp_test_df[
        gene_snp_test_df['gene_symbol'] == gene_ids.iloc[idx]
    ]
    # Subset genotype file for the significant SNPs
    genotype_df = get_genotype_df(residual_df, gene_snp_test_df)

    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene_symbol = df.gene_symbol
        gene_id = df.gene_id
        a1 = df.a1
        a2 = df.a2
        snp = df.snp_id
        gt = genotype_df[genotype_df.snp_id == snp][['sampleid', 'n_alt_alleles']]

        res_val = residual_df[['sampleid', gene_symbol]]
        test_df = res_val.merge(gt, on='sampleid', how='right')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        # Remove SNPs with NA values. This is due to some individuals having NA genotype values when
        # the residual df was created using the top eSNP as the conditioning SNP.
        test_df = test_df.dropna(axis=0, how='any')
        spearmans_rho, p = spearmanr(
            test_df['SNP'], test_df['residual'], nan_policy='omit'
        )
        return (gene_symbol, gene_id, a1, a2, snp, spearmans_rho, p)

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
    adjusted_spearman_df.to_parquet(output_path)

    return output_path


# Create click command line to enter dependency files
@click.command()
@click.option(
    '--significant-snps',
    required=True,
    help='A TSV of SNPs (as rows), \
        and significance levels, spearmans rho, and gene information as columns.',
)
@click.option(
    '--residuals',
    required=True,
    help='A CSV of gene residuals, with genes \
    as columns and samples as rows.',
)
@click.option(
    '--iterations', type=int, default=4, help='Number of iterations to perform'
)
@click.option(
    '--test-subset-genes',  # pylint: disable=too-many-locals
    type=int,
    help='Test with {test_subset_genes} genes, often = 5.',
)
def from_cli(
    significant_snps: str,
    residuals,
    iterations=4,
    test_subset_genes=None,
):

    # get chromosome and cell type info
    residual_list = residuals.split('/')
    residual_list = residual_list[:-1]
    chrom = residual_list[-1]
    cell_type = residual_list[-2]

    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name=f'Conditional analysis: {cell_type}/{chrom}',
        backend=backend,
        default_python_image=get_config()['workflow']['driver_image'],
    )

    main(
        batch=batch,
        job_prefix='',
        significant_snps_directory=significant_snps,
        residuals_path=residuals,
        chromosome=chrom,
        iterations=iterations,
        test_subset_genes=test_subset_genes,
        cell_type=cell_type,
        output_prefix=output_path(f'{cell_type}/{chrom}'),
    )


def main(
    batch: hb.Batch,
    *,
    job_prefix: str,
    cell_type: str,
    chromosome: str,
    significant_snps_directory: str,
    residuals_path,
    iterations: int = 4,
    output_prefix: str,
    test_subset_genes: int = 5,
    force: bool = False,
):
    """
    Creates a Hail Batch pipeline for calculating eQTLs using {iterations} iterations,
    scattered across the number of genes. Note, iteration 1 is completed in `generate_eqtl_spearan.py`.
    """

    # TODO: This isn't possible because residuals_path and significant_snps_directory_path
    #   are not available at runtime (they are intermediate outputs). This is currently not
    #   possible in Hail Batch. ** IF: this is the same number of scatters for this
    #   gene / chr pair as eqtl_spearman, we can calculate it up front in launch, and pass
    #   it down to here. If not, :
    #       - we need to work out another mechanism for determining it up-front
    #       - OR this must be a completely separate batch.
    #       - OR we wait for open batches
    n_genes = test_subset_genes or get_number_of_scatters(
        residuals_path=residuals_path,
        significant_snps_path=significant_snps_directory,
    )

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

        calc_resid_df_job = batch.new_python_job(
            f'{job_prefix}calculate-resid-df-iter-{iteration}'
        )
        calc_resid_df_job.depends_on(*sig_snps_dependencies)
        calc_resid_df_job.cpu(2)
        calc_resid_df_job.memory('8Gi')
        calc_resid_df_job.storage('2Gi')
        copy_common_env(calc_resid_df_job)

        round_dir = os.path.join(output_prefix, f'round{iteration}/')

        # this calculate_residual_df will also write an output that is useful later
        # if we assign previous_residual_path to the result of the call, we get the job dependency we actually want
        # as opposed to saying the following, which means there's actually no job dependency (bad):
        #    previous_residual_path = <output_path>
        previous_residual_path = calc_resid_df_job.call(
            calculate_residual_df,
            previous_residual_path,
            previous_sig_snps_directory,
            output_path=os.path.join(round_dir, f'residual_results.csv'),
            force=force,
        )

        sig_snps_component_paths = []

        new_sig_snps_directory = os.path.join(round_dir, 'sigsnps/')
        sink_jobs = []
        for gene_idx in range(n_genes):
            j = batch.new_python_job(
                name=f'{job_prefix}calculate_spearman_iter_{iteration}_job_{gene_idx}'
            )
            j.cpu(2)
            j.memory('8Gi')
            j.storage('2Gi')
            j.image(MULTIPY_IMAGE)
            j.depends_on(*sig_snps_dependencies)
            copy_common_env(j)
            gene_result_path: hb.resource.PythonResult = j.call(
                run_computation_in_scatter,
                iteration,
                gene_idx,
                previous_residual_path,
                previous_sig_snps_directory,
                cell_type,
                round_outputdir=round_dir,
                sigsnps_outputdir=new_sig_snps_directory,
                force=force,
            )
            sig_snps_component_paths.append(gene_result_path)
            sink_jobs.append(j)

        previous_sig_snps_directory = new_sig_snps_directory
        sig_snps_dependencies = sink_jobs

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


def fake_dependency(argument):
    return argument


if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    from_cli()
