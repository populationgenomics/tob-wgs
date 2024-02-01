"""Read in correlation results and claculate ld values"""

import click
import pandas as pd

import hail as hl


@click.command()
@click.option('--input-path', help='A path prefix of where to input files are located')
def query(input_path: str):
    """calculate LD values"""

    hl.init(default_reference='GRCh38')

    tmp_dir = input_path.replace(
        input_path.split('/')[2],
        input_path.split('/')[2] + '-tmp',
    )
    mt_path = f'{tmp_dir}/genotype_table.mt/'
    mt = hl.read_matrix_table(mt_path)
    mt = mt.annotate_rows(
        global_bp=hl.locus(mt.locus.contig, mt.locus.position).global_position(),
    )
    # if running on entire genome, concatenate all significant
    # correlation matrices together, then read in as one df
    # example: https://stackoverflow.com/questions/20906474/import-multiple-csv-files-into-pandas-and-concatenate-into-one-dataframe
    significant_snps_path = f'{input_path}correlation_results.csv'
    significant_snps_df = pd.read_csv(
        significant_snps_path,
        sep=' ',
        skipinitialspace=True,
    )
    t = hl.Table.from_pandas(significant_snps_df)
    # only keep rows whose FDR is < 0.05
    t = t.filter(t.fdr < 0.05)  # noqa: PLR2004
    print(f'Printing table: {t.show()}')
    t = t.key_by('global_bp')
    # filter mt to positions which are in significant_snps table
    mt = mt.filter_rows(hl.is_defined(t[mt.global_bp]))
    # add row index to be able to remap
    mt = mt.add_row_index()
    # turn matrix into table and save, in order to reference row idx
    mt_path = f'{input_path}significant_snps.mt'
    mt.write(mt_path)
    # perform ld calculation between all pairs of variants
    # within two megabases (1 megabase on either side)
    ld = hl.ld_matrix(mt.GT.n_alt_alleles(), mt.locus, radius=2e6)
    # only calculate the upper triangle
    ld = ld.sparsify_triangle()
    table = ld.entries()
    # replace row idx with global_bp
    table = table.rename({'i': 'row_idx'}).key_by('row_idx')
    mt = mt.key_rows_by('row_idx')
    table = (
        table.annotate(i=mt.rows()[table.row_idx].global_bp).key_by().drop('row_idx')
    )
    table = table.rename({'j': 'row_idx'}).key_by('row_idx')
    table = (
        table.annotate(j=mt.rows()[table.row_idx].global_bp).key_by().drop('row_idx')
    )
    # export as pandas table and save as csv
    table = table.to_pandas()
    # save table
    ld_filename = f'{input_path}ld_matrix.ht'
    # ld_filename = output_path(f'ld_matrix.csv', 'analysis')   # noqa: ERA001
    table.to_csv(ld_filename, index=False)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
