"""
Copy a sample-QC metadata TSV file generated by the joint-calling pipeline into the test bucket.
The idea is that for interactive exploration of QC summary metrics and design of plots, it's safe
to allow using the full set of samples, because the metadata doesn't contain any actual genomics
variantion data. 

Note that the test bucket also contains gs://cpg-tob-wgs-test/joint-calling/v2/meta.tsv, which
is symmetrical to gs://cpg-tob-wgs-main/joint-calling/v2/meta.tsv and limited to samples in the
test subset.
"""

import os
import subprocess
import click

output = os.getenv('OUTPUT')
assert output and output.startswith('gs://cpg-tob-wgs-test/')


@click.command()
@click.option('--version', help='Joint-calling run version', required=True)
def main(version):
    """Entry point"""
    subprocess.run(
        [
            'gsutil',
            'cp',
            f'gs://cpg-tob-wgs-main/joint-calling/{version}/meta.tsv',
            f'gs://cpg-tob-wgs-test/metadata/joint-calling/{version}/meta.tsv',
        ],
        check=False,
    )


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
