#!/usr/bin/env python3

"""
Use VEP using a dataproc cluster.
Taken from Matt Welland's script, run_vep_help.py

"""


import click
import hailtop.batch as hb
from analysis_runner import dataproc
from cpg_utils.hail_batch import get_config, remote_tmpdir

@click.command()
@click.option('--script', 'script', help='path to VEP main script')
@click.option('--mt', required=True, help='Hail matrix table to run VEP on')
def main(script: str, mt: str):
    """
    runs a script inside dataproc to execute VEP
    :param script: str, the path to the VEP main script
    """
    
    config = get_config()
    service_backend = hb.ServiceBackend(
        billing_project=config['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )

    # create a hail batch
    batch = hb.Batch(name='run_vep_in_dataproc_cluster', backend=service_backend)

    dataproc.hail_dataproc_job(
        batch=batch,
        worker_machine_type='n1-highmem-8',
        script=f'{script} --mt {mt}',
        max_age='12h',
        init=[
            'gs://cpg-reference/hail_dataproc/install_common.sh',
            'gs://cpg-reference/vep/vep-GRCh38.sh',
        ],
        job_name='run_vep',
        num_secondary_workers=20,
        num_workers=2,
        cluster_name='run vep',
    )

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
