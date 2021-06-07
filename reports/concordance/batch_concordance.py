#!/usr/bin/env python3

"""Concordance between SNPchip and WGS samples"""

import os
import hailtop.batch as hb


def concordance(batch, snpmt, wgsmt):
    """
    Concordance between SNPchip and WGS samples
    """
    conc = batch.new_job(name='run-concordance')
    conc.image('pdiakumis/concordance:0.1.2')
    conc.command(
        f"""
        set -e
        concordance \
          --snp {snpmt} \
          --wgs {wgsmt} \
          --res {conc.res} \
          --html {conc.html}
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
    snp = f'{BUCKET}/snpchip/v1/snpchip_grch38.mt'
    wgs = f'{BUCKET}/mt/test-v1-raw.mt'
    HTML = 'concordance_snpchip_with_wgs_chr22.html'
    RES = 'concordance_snpchip_with_wgs_chr22.tsv'
    concordance = concordance(b, snp, wgs)
    b.write_output(concordance.html, f'{BUCKET}-web/concordance/v1/{HTML}')
    b.write_output(concordance.res, f'{BUCKET}-web/concordance/v1/{RES}')
    b.run(dry_run=False)
