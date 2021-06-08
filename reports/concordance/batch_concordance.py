#!/usr/bin/env python3

"""Concordance between SNPchip and WGS samples"""

import os
import hailtop.batch as hb


def concordance(batch, snpmt, wgsmt, cpu):
    """
    Concordance between SNPchip and WGS samples
    """
    conc = batch.new_job(name='run-concordance')
    conc.image('pdiakumis/concordance:0.1.3')
    conc.cpu(8)
    conc.memory('16G')
    conc.storage('100G')
    conc.command(
        f"""
        set -e
        concordance \
          --snp {snpmt} \
          --wgs {wgsmt} \
          --res {conc.ores} \
          --html {conc.ohtml} \
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
    CPU = 8
    HTML = 'concordance_snpchip_with_wgs_chr22.html'
    RES = 'concordance_snpchip_with_wgs_chr22.tsv'
    concordance = concordance(b, SNP, WGS, CPU)
    b.write_output(concordance.ohtml, f'{BUCKET}-web/concordance/v1/{HTML}')
    b.write_output(concordance.ores, f'{BUCKET}-web/concordance/v1/{RES}')
    b.run(dry_run=False)
    service_backend.close()
