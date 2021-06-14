#!/usr/bin/env python3

"""Concordance between SNPchip and WGS samples"""

import os
import hailtop.batch as hb


def concordance(batch, snpmt, wgsmt, cpu):
    """
    Concordance between SNPchip and WGS samples
    """
    conc = batch.new_job(name='run-concordance')
    conc.image('pdiakumis/concordance:0.1.10')
    conc.cpu(cpu)
    conc.memory('highmem')
    conc.storage('100G')
    conc.command(
        f"""
        set -e
        concordance \
          --snp {snpmt} \
          --wgs {wgsmt} \
          --res_samples {conc.res_samples_tsv} \
          --res_sites {conc.res_sites_tsv} \
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

    BUCKET = 'gs://cpg-tob-wgs-test'
    SNP = f'{BUCKET}/snpchip/v1/snpchip_grch38.mt'
    WGS = f'{BUCKET}/mt/test-v1-raw.mt'
    CPU = 16
    PREFIX = 'concordance_TOB1524_chr22'
    HTML = f'{PREFIX}.html'
    concordance = concordance(b, SNP, WGS, CPU)
    b.write_output(concordance.html, f'{BUCKET}-web/concordance/v1/{HTML}')
    b.write_output(
        concordance.res_samples_tsv, f'{BUCKET}/concordance/v1/{PREFIX}_samples.tsv'
    )
    b.write_output(
        concordance.res_sites_tsv, f'{BUCKET}/concordance/v1/{PREFIX}_sites.tsv'
    )
    b.run(dry_run=False)
    service_backend.close()
