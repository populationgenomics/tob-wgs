"""
Launch analysis runner for all cell types and chromosomes

For example:

    python3 scripts/hail_batch/eqtl_hail_batch/launch_generate_eqtl_spearman.py \
        --input-path "gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files" \
        --dry-run \
        --chromosomes '1 2'
"""

import logging
import os
import click
import hail as hl
import hailtop.batch as hb
from cpg_utils.hail_batch import copy_common_env, remote_tmpdir
from cpg_utils.git import (
    get_git_commit_ref_of_current_repository,
    get_organisation_name_from_current_directory,
    get_repo_name_from_current_directory,
    prepare_git_job,
)
from cpg_utils.config import get_config
from google.cloud import storage


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
    help=('A path prefix of where input files are located. eg: gs://MyBucket/folder. '),
)
@click.option('--dry-run', is_flag=True, help='Just check if files exist')
def submit_eqtl_jobs(chromosomes, input_path, dry_run=False, cell_types=None):
    """Run association script for all chromosomes and cell types"""

    chromosomes = chromosomes.split(' ')
    assert isinstance(chromosomes, list)

    cs_client = storage.Client()

    bucket_name = input_path.split('gs://')[1].split('/')[0]
    bucket = cs_client.get_bucket(bucket_name)
    bucket_path = input_path.split(f'gs://{bucket_name}/')[-1]

    if cell_types is None or len(cell_types) == 0:
        # not provided (ie: use all cell types)
        # we can infer the cell types from the 'expression_files'
        # subdirectory of the input_path
        # eg: {cell_type}_expression.tsv

        path_to_expression_files = os.path.join(bucket_path, 'expression_files')
        logging.info(f'Going to fetch cell types from {path_to_expression_files}')
        blobs = bucket.list_blobs(prefix=path_to_expression_files + '/', delimiter='/')
        ending = '_expression.tsv'

        cell_types = [
            os.path.basename(b.name)[: -len(ending)]
            for b in blobs
            if b.name.endswith('_expression.tsv')
        ]
        logging.info(f'Found {len(cell_types)} cell types: {cell_types}')

    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(name='eqtl_spearman', backend=backend)

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
                lambda x: not hl.hadoop_exists(x), files_to_check
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
                    lambda x: not hl.hadoop_exists(x), files_to_check
                )
                for file in files_that_are_missing:
                    logging.error(f'File {file} is missing')
            else:

                job = batch.new_job(f'{cell_type}-chr{chromosome}')
                copy_common_env(job)
                job.image(get_config()['workflow']['driver_image'])

                # check out a git repository at the current commit
                prepare_git_job(
                    job=job,
                    organisation=get_organisation_name_from_current_directory(),
                    repo_name=get_repo_name_from_current_directory(),
                    commit=get_git_commit_ref_of_current_repository(),
                )

                keys = os.path.join(input_path, 'OneK1K_CPG_IDs.tsv')
                job.command(
                    f'python3 scripts/eqtl_hail_batch/generate_eqtl_spearman.py '
                    f'--expression {expression} '
                    f'--covariates {covariates} '
                    f'--geneloc {geneloc} '
                    f'--keys {keys} '
                )

    batch.run(wait=False)


if __name__ == '__main__':
    submit_eqtl_jobs()  # pylint: disable=no-value-for-parameter
