"""Concordance between SNPchip and WGS samples"""

import hailtop.batch as hb


def concordance(b, snpmt, wgsmt):
    """
    Concordance between SNPchip and WGS samples
    """
    c = b.new_job(name='run-concordance')
    c.image('pdiakumis/concordance:0.1')
    c.command(
        f"""
        set -e
        concordance \
          --snp {snpmt} \
          --wgs {wgsmt} \
          --html {c.ofile}
        """
    )
    return c


if __name__ == '__main__':
    batch = hb.Batch(name='concordance')
    BUCKET = 'gs://cpg-tob-wgs-test'
    snp = batch.read_input(f'{BUCKET}/snpchip/v1/snpchip_grch38.mt')
    wgs = batch.read_input(f'{BUCKET}/mt/test-v1-raw.mt')
    html = f'concordance_snpchip_with_wgs.html'
    concordance = concordance(batch, snp, wgs)
    batch.write_output(concordance.ofile, f'{BUCKET}-web/concordance/v1/{html}')
    batch.run(dry_run=False)
