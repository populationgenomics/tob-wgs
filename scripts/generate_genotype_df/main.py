"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from cpg_utils.hail import remote_tmpdir
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('CPG_DATASET'), remote_tmpdir=remote_tmpdir()
)

batch = hb.Batch(name='export_genotype_data', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'generate_genotype_df.py',
    max_age='4h',
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='export_genotype_data',
)

batch.run()
