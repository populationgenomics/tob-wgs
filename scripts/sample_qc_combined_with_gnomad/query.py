"""Generates loadings, scores, and eigenvalues for the HGDP,1KG, and tob-wgs dataset"""

import click
import hail as hl


@click.command()
@click.option('--test', is_flag=True)
@click.option('--output', help='GCS output path', required=True)
def query(test: bool, output: str):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hgdp1kg_tobwgs_joined_path = 'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v0/hgdp1kg_tobwgs_joined_all_samples.mt/'
    if test:
        hgdp1kg_tobwgs_joined_path = hgdp1kg_tobwgs_joined_path.replace(
            'gs://cpg-tob-wgs-analysis', 'gs://cpg-tob-wgs-test'
        )

    hgdp1kg_tobwgs_joined = hl.read_matrix_table(hgdp1kg_tobwgs_joined_path)
    qc_ht = hl.sample_qc(hgdp1kg_tobwgs_joined, name='sample_qc')
    qc_ht.write(output)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
