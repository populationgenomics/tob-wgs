""" Copy VEP reference data to a local bucket
"""

import hail as hl
import hailtop.batch as hb

hl.init(default_reference='GRCh38')

service_backend = hb.ServiceBackend(
    billing_project='vladislavsavelyev-trial', bucket='gs://playground-au/hail'
)

b = hb.Batch(name='dataproc example', backend=service_backend)
j = b.new_job('Copy VEP data')
j.command('gcloud -q auth activate-service-account --key-file=/gsa-key/key.json')
j.command('gcloud config set project hail-295901')
j.command(
    'gsutil -q -u hail-295901 cp -r gs://gnomad-public-requester-pays/resources/context/grch38_context_vep_annotated.ht gs://cpg-reference/hg38/v0/grch38_context_vep_annotated.ht'
)
j.image(
    'australia-southeast1-docker.pkg.dev/analysis-runner/images/dataproc:hail-0.2.64'
)
b.run()
