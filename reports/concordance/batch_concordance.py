#!/usr/bin/env python3

"""Concordance between SNPchip and WGS samples"""

import os
import hailtop.batch as hb

CONCORDANCE_IMG = (
    'australia-southeast1-docker.pkg.dev/peter-dev-302805/test/concordance:0.1.16'
)


def concordance(batch, snpmt, wgsmt, samples, chrom, cpu):
    """
    Concordance between SNPchip and WGS samples
    """
    conc = batch.new_job(name='run-concordance')
    conc.image(CONCORDANCE_IMG)
    conc.cpu(cpu)
    conc.memory('lowmem')
    conc.storage('100G')
    input_samples = batch.read_input(samples)
    conc.command(
        f"""
        set -e
        concordance \
          --snp {snpmt} \
          --wgs {wgsmt} \
          --samples {input_samples} \
          --chrom {chrom} \
          --res_samples {conc.res_samples_tsv} \
          --html {conc.html} \
          --cpu {cpu}
        """
    )
    return conc


if __name__ == '__main__':
    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )
    b = hb.Batch(backend=service_backend, name='concordance')

    BUCKET = 'gs://cpg-tob-wgs-main'
    SNP = f'{BUCKET}/snpchip/v1/snpchip_grch38.mt'
    WGS = f'{BUCKET}/mt/v5.1-nonref.mt'
    SAMPLES = 'gs://cpg-tob-wgs-test/pdiakumis/concordance/samples_to_keep.tsv'
    CHROM = 'chr22'
    CPU = 32
    PREFIX = f'v5_subset_samples_{CHROM}'
    HTML = f'{PREFIX}.html'
    concordance = concordance(b, SNP, WGS, SAMPLES, CHROM, CPU)
    b.write_output(concordance.html, f'{BUCKET}-web/concordance/v1/{HTML}')
    b.write_output(
        concordance.res_samples_tsv, f'{BUCKET}/concordance/v1/{PREFIX}_samples.tsv'
    )
    b.run(dry_run=False)
    service_backend.close()
