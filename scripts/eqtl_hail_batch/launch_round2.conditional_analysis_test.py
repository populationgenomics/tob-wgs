"""
Launch analysis runner for all cell types and chromosomes

For example:

    python3 scripts/hail_batch/eqtl_hail_batch/launch_round2.conditional_analysis_test.py \
        --input-path "gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files" \
        --dry-run \
        --output-dir gs://cpg-tob-wgs-test/eqtl_output \
        --chromosomes '1 2'
"""

import logging
import os
import click

from analysis_runner.cli_analysisrunner import run_analysis_runner
from google.cloud import storage


# def run_analysis_runner(*args, **kwargs):
#     """Run analysis runner for all cell types and chromosomes"""
#     print(args, kwargs)


@click.command()
@click.option(
    '--cell-types',
    default=None,
    help='List of cell types to test',
)
@click.option(
    '--chromosomes',
    required=True,
    help='List of chromosome numbers to run eQTL analysis on. Space separated, as one argument',  # noqa: E501; pylint: disable=line-too-long
)
@click.option(
    '--input-path',
    required=True,
    help=(
        'A path prefix of where input files are located, eg: gs://MyBucket/folder/. '
        'We expect 5 folders to exist in this subdirectory: '
        'expression_files, gene_location_files, genotype_files, snp_location_files, covariates_files'  # noqa: E501; pylint: disable=line-too-long
    ),
)
@click.option(
    '--first-round-path',
    required=True,
    help=(
        'A path prefix of where the output from the first round of esnp '
        'results are located, eg: gs://MyBucket/folder/. '
    ),
)
@click.option(
    '--output-dir',
    required=True,
    help='A path of where to output files, eg: gs://MyBucket/output-folder/',
)
@click.option(
    '--test-subset-genes',  # pylint: disable=too-many-locals
    type=int,
    help='Test with {test_subset_genes} genes, often = 5.',
)
@click.option('--dry-run', is_flag=True, help='Just check if files exist')
def submit_eqtl_jobs(
    chromosomes,
    input_path,
    first_round_path,
    output_dir,
    test_subset_genes=None,
    dry_run=False,
    cell_types=None,
):
    """Run association script for all chromosomes and cell types"""

    assert output_dir.startswith('gs://') and input_path.startswith('gs://')
    chromosomes = chromosomes.split(' ')
    assert isinstance(chromosomes, list)

    cs_client = storage.Client()

    bucket_name = input_path.split('gs://')[1].split('/')[0]
    bucket = cs_client.get_bucket(bucket_name)
    bucket_path = input_path.split(f'gs://{bucket_name}/')[-1]

    def file_exists(files):
        file_status = storage.Blob(
            bucket=bucket, name=files.split(bucket_name)[-1].strip('/')
        ).exists(cs_client)
        return file_status

    if cell_types is None:
        # not provided (ie: use all cell types)
        # we can infer the cell types from the 'expression_files'
        # subdirectory of the input_path
        # eg: {cell_type}_expression.tsv

        path_to_expression_files = os.path.join(bucket_path, 'expression_files')
        logging.info(f'Going to fetch cell types from {path_to_expression_files}')
        blobs = bucket.list_blobs(prefix=path_to_expression_files + '/', delimiter='/')

        cell_types = [
            os.path.basename(b.name)[:-15]
            for b in blobs
            if b.name.endswith('_expression.tsv')
        ]
        logging.info(f'Found {len(cell_types)} cell types: {cell_types}')

    for cell_type in cell_types:
        expression = os.path.join(
            input_path, 'expression_files', f'{cell_type}_expression.tsv'
        )
        covariates = os.path.join(
            input_path, 'covariates_files', f'{cell_type}_peer_factors.tsv'
        )

        if dry_run:
            files_to_check = [expression, covariates]
            files_that_are_missing = filter(
                lambda x: not file_exists(x), files_to_check
            )
            for file in files_that_are_missing:
                logging.error(f'File {file} is missing')

        for chromosome in chromosomes:
            residuals = os.path.join(
                first_round_path, cell_type, f'chr{chromosome}', f'log_residuals.tsv'
            )
            significant_snps = os.path.join(
                first_round_path,
                cell_type,
                f'chr{chromosome}',
                f'correlation_results.csv',
            )
            genotype = os.path.join(
                input_path, 'genotype_files', f'tob_genotype_chr{chromosome}.tsv'
            )

            if dry_run:
                # check all files exist before running
                files_to_check = [genotype, residuals]
                files_that_are_missing = filter(
                    lambda x: not file_exists(x), files_to_check
                )
                for file in files_that_are_missing:
                    logging.error(f'File {file} is missing')
            else:

                output_prefix = os.path.join(
                    output_dir, f'eqtl_results_{cell_type}', f'chr{chromosome}'
                )
                # add in sampleid keys
                keys = os.path.join(input_path, 'OneK1K_CPG_IDs.tsv')
                # The analysis-runner output path doesn't want the BUCKET specified,
                # so let's remove it from the output_prefix
                analysis_runner_output_path = output_prefix[5:].partition('/')[-1]
                run_analysis_runner(
                    description=f'eqtl_spearman_{cell_type}_chr{chromosome}',
                    dataset='tob-wgs',
                    access_level='test',
                    output_dir=analysis_runner_output_path,
                    # commit, sha and cwd can be inferred automatically
                    script=[
                        'round2.conditional_analysis_test.py',
                        *('--residuals', residuals),
                        *('--significant-snps', significant_snps),
                        *('--genotype', genotype),
                        *('--output-prefix', output_prefix),
                        *(
                            ['--test-subset-genes', str(test_subset_genes)]
                            if test_subset_genes
                            else []
                        ),
                        *('--keys', keys),
                    ],
                )


if __name__ == '__main__':
    submit_eqtl_jobs()  # pylint: disable=no-value-for-parameter
