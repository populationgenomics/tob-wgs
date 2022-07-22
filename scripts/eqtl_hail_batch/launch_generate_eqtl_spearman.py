"""
Launch analysis runner for all cell types and chromosomes

For example:

    python3 scripts/hail_batch/eqtl_hail_batch/launch_generate_eqtl_spearman.py \
        --input-path "gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files" \
        --dry-run \
        --output-dir gs://cpg-tob-wgs-test/eqtl_output \
        --chromosomes '1 2'
"""

import logging
import os
import click
import subprocess
import hailtop.batch as hb
# from analysis_runner.cli_analysisrunner import run_analysis_runner
from google.cloud import storage
from cpg_utils.git import (
  get_git_commit_ref_of_current_repository,
  get_organisation_name_from_current_directory,
  get_repo_name_from_current_directory,
  prepare_git_job,
)


@click.command()
@click.option(
    '--cell-types',
    default=None,
    multiple=True,
    help='List of cell types to test. All available cell types can \
        be found in \
        `gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/expression_files/`',
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
    ),
)
@click.option(
    '--output-dir',
    required=True,
    help='A path of where to output files, eg: gs://MyBucket/output-folder/',
)
@click.option('--dry-run', is_flag=True, help='Just check if files exist')
def submit_eqtl_jobs(
    chromosomes, input_path, output_dir, dry_run=False, cell_types=None
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
        ending = '_expression.tsv'

        cell_types = [
            os.path.basename(b.name)[:-len(ending)]
            for b in blobs
            if b.name.endswith('_expression.tsv')
        ]
        logging.info(f'Found {len(cell_types)} cell types: {cell_types}')

    for cell_type in cell_types:
        expression = os.path.join(
            input_path, 'expression_files', f'{cell_type}_expression.tsv'
        )
        covariates = os.path.join(
            input_path, 'covariates_files', f'{cell_type}_peer_factors_file.txt'
        )

        if dry_run:
            files_to_check = [expression, covariates]
            files_that_are_missing = filter(
                lambda x: not file_exists(x), files_to_check
            )
            for file in files_that_are_missing:
                logging.error(f'File {file} is missing')

        for chromosome in chromosomes:
            geneloc = os.path.join(
                input_path, 'gene_location_files', f'GRCh38_geneloc_chr{chromosome}.tsv'
            )

            if dry_run:
                # check all files exist before running
                files_to_check = [geneloc]
                files_that_are_missing = filter(
                    lambda x: not file_exists(x), files_to_check
                )
                for file in files_that_are_missing:
                    logging.error(f'File {file} is missing')
            else:

                batch = hb.Batch('run-eqtl-analysis')
                job = batch.new_job('checkout_repo')

                # check out a git repository at the current commit
                prepare_git_job(
                    job=job,
                    organisation=get_organisation_name_from_current_directory(),
                    repo_name=get_repo_name_from_current_directory(),
                    commit=get_git_commit_ref_of_current_repository(),
                )

                keys = os.path.join(input_path, 'OneK1K_CPG_IDs.tsv')
                output_prefix = os.path.join(
                    output_dir, f'{cell_type}', f'chr{chromosome}'
                )
                subprocess.check_output(
                    ['python3',
                     'generate_eqtl_spearman.py',
                     *('--expression', expression),
                     *('--covariates', covariates),
                     *('--geneloc', geneloc),
                     *('--keys', keys),
                     *('--output-prefix', output_prefix),
                     ]
                     )


if __name__ == '__main__':
    submit_eqtl_jobs()  # pylint: disable=no-value-for-parameter
