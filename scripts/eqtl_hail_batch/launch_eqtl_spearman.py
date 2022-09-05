"""
Create a Hail Batch workflow for all the EQTL analysis, including:

- generate_eqtl_spearman
- conditional_analysis rounds

For example:

    python3 scripts/hail_batch/eqtl_hail_batch/launch_eqtl_spearman.py \
        --input-files-prefix "gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files" \
        --chromosomes '1 2' \
        --genes B_intermediate
"""

import logging
import os
from collections import defaultdict

import click
import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    remote_tmpdir,
    dataset_path,
    output_path,
    reference_path,
    copy_common_env,
)
from google.cloud import storage

from generate_eqtl_spearman import main as generate_eqtl_spearman
from genotype_info import filter_joint_call_mt

DEFAULT_JOINT_CALL_TABLE_PATH = (dataset_path('mt/v7.mt/'),)
DEFAULT_FREQUENCY_TABLE_PATH = (
    dataset_path('joint-calling/v7/variant_qc/frequencies.ht/', 'analysis'),
)
DEFAULT_VEP_ANNOTATION_TABLE_PATH = (
    dataset_path('tob_wgs_vep/104/vep104.3_GRCh38.ht/'),
)
DEFAULT_GENCODE_GTF_PATH = 'gs://cpg-reference/gencode/gencode.v38.annotation.gtf.bgz'  # reference_path('gencode_gtf'),


@click.command()
@click.option(
    '--input-files-prefix',
    required=True,
    help='A path prefix of where input files are located. eg: gs://MyBucket/folder. '
    'If a relative path is given, it will be from the output-path',
)
@click.option(
    '--cell-types',
    default=None,
    multiple=True,
    help='List of cell types to test. All available cell types can be found in '
    '`gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/expression_files/`',
)
@click.option(
    '--chromosomes',
    required=True,
    help='List of chromosome numbers to run eQTL analysis on. '
    'Space separated, as one argument',
)
@click.option('--force', is_flag=True, help='Skip checkpoints')
@click.option('--dry-run', is_flag=True, help='Just check if files exist')
def from_cli(chromosomes, input_files_prefix, cell_types):
    chromosomes_list = chromosomes.split(' ')

    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(name='eqtl_spearman', backend=backend)

    _ = main(
        batch,
        input_files_prefix=input_files_prefix,
        chromosomes=chromosomes_list,
        cell_types=cell_types,
    )
    batch.run(wait=False, dry_run=True)


def main(
    batch: hb.Batch,
    *,
    input_files_prefix: str,
    chromosomes: list[str],
    cell_types: list[str] = None,
    joint_call_table_path: str = DEFAULT_JOINT_CALL_TABLE_PATH,
    frequency_table_path: str = DEFAULT_FREQUENCY_TABLE_PATH,
    vep_annotation_table_path: str = DEFAULT_VEP_ANNOTATION_TABLE_PATH,
    gencode_gtf_path: str = DEFAULT_GENCODE_GTF_PATH,
    force: bool=False
):
    """Run association script for all chromosomes and cell types"""

    if not any(
        input_files_prefix.startswith(prefix) for prefix in ("gs://", '/', 'https://')
    ):
        input_files_prefix = output_path(input_files_prefix)

    if not (isinstance(chromosomes, list) and len(chromosomes) > 0):
        raise ValueError('Must specify at least 1 chromosome as a list')

    if cell_types is None or len(cell_types) == 0:
        # not provided (ie: use all cell types)
        cell_types = get_cell_types_from(dataset_path('expression_files'))
        logging.info(f'Found {len(cell_types)} cell types: {cell_types}')

    # ideally this would come from metamist :(
    keys_tsv_path = os.path.join(input_files_prefix, 'OneK1K_CPG_IDs.tsv')

    outputs = defaultdict(dict)

    # do the genotype_info stuff
    filter_mt_job = batch.new_python_job(f'filter_mt')
    copy_common_env(filter_mt_job)
    filtered_mt_path = filter_mt_job.call(
        filter_joint_call_mt,
        keys_path=keys_tsv_path,
        joint_mt_path=joint_call_table_path,
        frequency_table_path=frequency_table_path,
        vep_annotation_path=vep_annotation_table_path,
        output_path=output_path('genotype_table.mt', 'tmp'),
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

            outputs[cell_type][chromosome] = generate_eqtl_spearman(
                batch=batch,
                # constants
                force=force,
                job_prefix=f'{cell_type}_chr{chromosome}_',
                cell_type=cell_type,
                chromosome=chromosome,
                output_prefix=output_path(f'{cell_type}/{chromosome}'),
                # derived paths
                filtered_mt_path=filtered_mt_path,
                gencode_gtf_path=gencode_gtf_path,
                expression_tsv_path=expression_tsv_path,
                covariates_tsv_path=covariates_tsv_path,
                geneloc_tsv_path=os.path.join(
                    input_files_prefix,
                    'gene_location_files',
                    f'GRCh38_geneloc_chr{chromosome}.tsv',
                ),
            )


def get_cell_types_from(expression_files_dir: str):
    """
    we can infer the cell types from the 'expression_files' subdirectory
        eg: {cell_type}_expression.tsv
    """
    logging.info(f'Going to fetch cell types from {expression_files_dir}')

    ending = '_expression.tsv'

    if expression_files_dir.startswith("gs://"):
        cs_client = storage.Client()
        bucket_name, bucket_path = expression_files_dir.split('gs://')[1].split(
            '/', maxsplit=1
        )
        bucket = cs_client.get_bucket(bucket_name)

        blobs = bucket.list_blobs(prefix=bucket_path + '/', delimiter='/')
        return [
            os.path.basename(b.name)[: -len(ending)]
            for b in blobs
            if b.name.endswith(ending)
        ]

    if os.path.exists(expression_files_dir):
        return [
            os.path.basename(fn)[: -len(ending)]
            for fn in os.listdir(expression_files_dir)
            if fn.endswith(ending)
        ]

    raise ValueError(f'Unrecognised expression_files_dir type: {expression_files_dir}')


if __name__ == '__main__':

    backend = None
    # backend = hb.ServiceBackend(
    #     billing_project=get_config()['hail']['billing_project'],
    #     remote_tmpdir=remote_tmpdir(),
    # )
    batch = hb.Batch(name='eqtl_spearman', backend=backend)

    _ = main(
        batch,
        input_files_prefix='gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files',
        chromosomes=['22'],
        cell_types=['B_intermediate'],
    )
    batch.run(wait=False, dry_run=True)
