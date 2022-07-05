"""Read in correlation results and claculate ld values"""

import hail as hl
import click
import pandas as pd


@click.command()
@click.option('--input-path', help='A path prefix of where to input files are located')
def query(input_path):
    """calculate LD values"""

    hl.init(default_reference='GRCh38')

    tmp_dir = input_path.replace(input_path.split('/')[2], input_path.split('/')[2] + '-tmp')
    mt_path = f'{tmp_dir}/genotype_table.mt/'
    mt = hl.read_matrix_table(mt_path)
    mt = mt.annotate_rows(
            global_bp=hl.locus(
                mt.locus.contig, mt.locus.position
            ).global_position(),
        )
    significant_snps_path = f'{input_path}correlation_results.csv'
    significant_snps_df = pd.read_csv(
        significant_snps_path, sep=' ', skipinitialspace=True
    )
    t = hl.Table.from_pandas(significant_snps_df)
    # only keep rows whose FDR is < 0.05
    # t = t.filter(t.fdr < 0.05)
    print(f'Printing table: {t.show()}')
    t = t.key_by('global_bp')
    # filter mt to positions which are in significant_snps table
    significant_snps = mt.filter_rows(hl.is_defined(t[mt.global_bp]))
    # add row index to be able to remap
    significant_snps = significant_snps.add_row_index()
    # turn matrix into table and save, in order to reference row idx
    significant_snps_path = f'{input_path}significant_snps.mt'
    significant_snps.write(significant_snps_path)
    # perform ld calculation
    ld = hl.ld_matrix(significant_snps.GT.n_alt_alleles(), significant_snps.locus, radius=2e6)
    table = ld.entries()
    # filter out entries with an LD score less than 0.2
    table = table.filter(table.entry > 0.2)
    # save table
    ld_filename = f'{input_path}ld_matrix.ht'
    table.write(ld_filename)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
