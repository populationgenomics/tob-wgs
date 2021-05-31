"""Run hgdp_1kg_tob_wgs_pca.py using the analysis runner."""

import os
import hail as hl
import hailtop.batch as hb
from analysis_runner import dataproc

OUTPUT = os.getenv('OUTPUT')
assert OUTPUT

hl.init(default_reference='GRCh38')

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='sample_qc combined hgdp1kg tobwgs', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'query.py --output={OUTPUT}',
    max_age='24h',
    num_workers=50,
    packages=['click'],
    job_name='hgdp1kg-tobwgs-sample-qc',
)

batch.run()
