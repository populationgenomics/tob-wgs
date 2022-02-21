#!/usr/bin/env python3

"""Untar and unzip expression files"""

import os
import click
import hailtop.batch as hb

from analysis_runner.constants import GCLOUD_ACTIVATE_AUTH

@click.command()
@click.option(
    '--file',
    required=True,
    help='File to unzip and untar',
)
def main(file: str):
    """
    Untar and unzip files
    """
    dataset = os.getenv('DATASET')
    access_level = os.getenv('ACCESS_LEVEL')
    backend = hb.ServiceBackend(
        billing_project=dataset, bucket=f'cpg-{dataset}-{access_level}'
    )
    b = hb.Batch(name='untar-files', backend=backend)
    tarfile = b.read_input(file)
    j = b.new_job('untar-files')
    j.image(os.getenv('DRIVER_IMAGE'))
    j.command(f'tar xvfz {tarfile}')
    j.command(GCLOUD_ACTIVATE_AUTH)
    j.command(f"""gsutil cp expression_220218/* gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/expression_files/""")
    j.command(f"""gsutil mv expression_220218/B_intermediate_expression.tsv expression_220218/B_memory_expression.tsv expression_220218/B_naive_expression.tsv gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/expression_files/""")
    b.run(wait=False)

if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    main()
