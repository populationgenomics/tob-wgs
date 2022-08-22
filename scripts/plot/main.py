#!/usr/bin/env python3

"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc


backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='plot_figure', backend=backend)


dataproc.hail_dataproc_job(
    batch,
    f'plot_data.py',
    max_age='1h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh'],
    job_name=f'plot_data',
)

batch.run(wait=False)
