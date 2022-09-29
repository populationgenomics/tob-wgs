#!/usr/bin/env python3

import logging
import subprocess
import sys
import re
import click
import pandas as pd
from cloudpathlib import AnyPath
from cpg_utils import to_path
from cpg_utils.hail_batch import output_path

subprocess.run(
    [
        sys.executable,
        '-m',
        'pip',
        'install',
        'limix==3.0.4',
        'pandas_plink==2.2.9',
        'xarray==0.20.2',
    ],
    check=True,
)

import xarray as xr  # pylint: disable=wrong-import-position, import-error
from pandas_plink import (
    read_plink1_bin,
)  # pylint: disable=wrong-import-position, import-error
from limix.qc import (
    quantile_gaussianize,
)  # pylint: disable=wrong-import-position, import-error

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@click.command()
@click.option('--cell-type', required=True)  # 'B_naive'
@click.option('--gene-name', required=True)  # 'VPREB3'
@click.option(
    '--sample-mapping-file', required=True
)  # 'scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
@click.option(
    '--genotype-file-bed', required=True
)  # 'v0/plink_files/vpreb3_rare_promoter.bed'
@click.option(
    '--genotype-file-bim', required=True
)  # 'v0/plink_files/vpreb3_rare_promoter.bim'
@click.option(
    '--genotype-file-fam', required=True
)  # 'v0/plink_files/vpreb3_rare_promoter.fam'
@click.option(
    '--phenotype-file', required=True
)  # 'scrna-seq/grch38_association_files/expression_files/B_naive_expression.tsv'
@click.option('--kinship-file', required=True)  # 'v0/skat/grm_wide.csv'
@click.option(
    '--output-folder', required=False, default=''
)  # by default current directory, where you are running your script from
def main(
    cell_type: str,
    gene_name: str,
    sample_mapping_file: str,
    genotype_file_bed: str,
    genotype_file_bim: str,
    genotype_file_fam: str,
    phenotype_file: str,
    kinship_file: str,
    output_folder: str,  # do I need this??
):

    expression_filename = AnyPath(output_path(f'{gene_name}_{cell_type}.csv'))
    genotype_filename = AnyPath(output_path(f'{gene_name}_rare_regulatory.csv'))
    kinship_filename = AnyPath(output_path('kinship_common_samples.csv'))

    # region PHENOTYPE_FILE

    phenotype = pd.read_csv(phenotype_file, sep='\t', index_col=0)

    phenotype = xr.DataArray(
        phenotype.values,
        dims=['sample', 'gene'],
        coords={'sample': phenotype.index.values, 'gene': phenotype.columns.values},
    )

    # endregion PHENOTYPE_FILE

    # region KINSHIP_FILE

    ## read in GRM (genotype relationship matrix; kinship matrix)
    kinship = pd.read_csv(kinship_file, index_col=0)
    kinship.index = kinship.index.astype('str')
    assert all(kinship.columns == kinship.index)  # symmetric matrix, donors x donors

    kinship = xr.DataArray(
        kinship.values,
        dims=['sample_0', 'sample_1'],
        coords={'sample_0': kinship.columns, 'sample_1': kinship.index},
    )
    kinship = kinship.sortby('sample_0').sortby('sample_1')

    # endregion KINSHIP_FILE

    # region GENOTYPE_FILE

    # read in genotype file (plink format)
    # bed
    with to_path(genotype_file_bed).open('rb') as handle:
        data = handle.readlines()
    with open('temp.bed', 'wb') as handle:
        handle.writelines(data)
    # bim
    with to_path(genotype_file_bim).open('rb') as handle:
        data = handle.readlines()
    with open('temp.bim', 'wb') as handle:
        handle.writelines(data)
    # fam
    with to_path(genotype_file_fam).open('rb') as handle:
        data = handle.readlines()
    with open('temp.fam', 'wb') as handle:
        handle.writelines(data)
    # read
    geno = read_plink1_bin('temp.bed')

    # endregion GENOTYPE_FILE

    # region SAMPLE_MAPPING_FILE

    # this file will map different IDs (and OneK1K ID to CPG ID)
    sample_mapping = pd.read_csv(sample_mapping_file, sep='\t')

    ## samples with expression data
    donors_onek1k = sample_mapping['OneK1K_ID'].unique()
    donors_onek1k.sort()
    donors_exprs = sorted(
        set(list(phenotype.sample.values)).intersection(donors_onek1k)
    )
    logging.info(f'Number of unique donors with expression data: {len(donors_exprs)}')

    ## samples with genotype data
    donors_cpg = sample_mapping['InternalID'].unique()
    donors_cpg.sort()
    donors_geno = sorted(set(list(geno.sample.values)).intersection(donors_cpg))
    logging.info(f'Number of unique donors with genotype data: {len(donors_geno)}')

    ## samples with both (can this be done in one step?)
    sample_mapping1 = sample_mapping.loc[sample_mapping['OneK1K_ID'].isin(donors_exprs)]
    sample_mapping_both = sample_mapping1.loc[
        sample_mapping1['InternalID'].isin(donors_geno)
    ]
    donors_e = sample_mapping_both['OneK1K_ID'].unique()
    donors_g = sample_mapping_both['InternalID'].unique()
    assert len(donors_e) == len(donors_g)

    logging.info(f'Number of unique common donors: {len(donors_g)}')

    ## samples in kinship
    donors_e_short = [re.sub('.*_', '', donor) for donor in donors_e]
    donors_k = sorted(set(list(kinship.sample_0.values)).intersection(donors_e_short))

    # endregion SAMPLE_MAPPING_FILE

    # region SUBSET_FILES

    ###############
    #### phenotype
    phenotype = phenotype.sel(sample=donors_e)

    # select gene
    y = phenotype.sel(gene=gene_name)
    y = quantile_gaussianize(y)

    del phenotype  # delete to free up memory

    # make data frame to save as csv
    y_df = pd.DataFrame(
        data=y.values.reshape(y.shape[0], 1), index=y.sample.values, columns=[gene_name]
    )

    ###############
    #### genotype
    geno = geno.sel(sample=donors_g)

    # make data frame to save as csv
    data = geno.values
    geno_df = pd.DataFrame(data, columns=geno.snp.values, index=geno.sample.values)
    geno_df = geno_df.dropna(axis=1)

    # delete large files to free up memory
    del G

    ###############
    #### kinship
    kinship = kinship.sel(sample_0=donors_k, sample_1=donors_k)
    assert all(kinship.sample_0 == donors_k)
    assert all(kinship.sample_1 == donors_k)

    # make data frame to save as csv
    kinship_df = pd.DataFrame(
        kinship.values, columns=kinship.sample_0, index=kinship.sample_1
    )

    del kinship  # delete kinship to free up memory

    # endregion SUBSET_FILES

    # region SAVE_FILES

    with expression_filename.open('w') as ef:
        y_df.to_csv(ef, index=False)

    with genotype_filename.open('w') as gf:
        geno_df.to_csv(gf, index=False)

    with kinship_filename.open('w') as kf:
        kinship_df.to_csv(kf, index=False)

    # endregion SAVE_FILES


if __name__ == '__main__':  # --exclude
    main()  
