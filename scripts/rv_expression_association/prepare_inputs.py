#!/usr/bin/env python3

import click
import logging
import re
import subprocess
import sys
import pandas as pd
import xarray as xr
from cloudpathlib import AnyPath
from cpg_utils import to_path
from cpg_utils.hail_batch import output_path

subprocess.run([sys.executable, '-m', 'pip', 'install', 'limix==3.0.4', 'pandas_plink==2.2.9'], check=True)

from pandas_plink import read_plink1_bin   # pylint: disable=wrong-import-position
from limix.qc import quantile_gaussianize  # pylint: disable=wrong-import-position

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@click.command()
@click.option('--cell-type', required=True, help='More info here')
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
    output_folder: str,
):

    expression_filename = AnyPath(output_path(f'{gene_name}_{cell_type}.csv'))
    genotype_filename = AnyPath(output_path(f'{gene_name}_rare_regulatory.csv'))
    kinship_filename = AnyPath(output_path('kinship_common_samples.csv'))

    ##################### lint also not happy with these, should I care?
    ### phenotype file ##
    #####################

    phenotype = pd.read_csv(phenotype_file, sep='\t', index_col=0)

    phenotype = xr.DataArray(
        phenotype.values,
        dims=['sample', 'gene'],
        coords={'sample': phenotype.index.values, 'gene': phenotype.columns.values},
    )

    ####################
    ### kinship file ###
    ####################

    ## read in GRM (genotype relationship matrix; kinship matrix)
    K = pd.read_csv(kinship_file, index_col=0)
    K.index = K.index.astype('str')
    assert all(K.columns == K.index)  # symmetric matrix, donors x donors

    K = xr.DataArray(
        K.values,
        dims=['sample_0', 'sample_1'],
        coords={'sample_0': K.columns, 'sample_1': K.index},
    )
    K = K.sortby('sample_0').sortby('sample_1')

    #####################
    ### genotype files ##
    #####################

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
    G = read_plink1_bin('temp.bed')

    #################################
    ### sample mapping file (SMF) ###
    #################################

    # this file will map different IDs (and OneK1K ID to CPG ID)
    sample_mapping = pd.read_csv(sample_mapping_file, sep='\t')

    ## samples with expression data
    donors_onek1k = sample_mapping['OneK1K_ID'].unique()
    donors_onek1k.sort()
    donors_exprs = sorted(
        set(list(phenotype.sample.values)).intersection(donors_onek1k)
    )
    logging.info(
        'Number of unique donors with expression data: {}'.format(len(donors_exprs))
    )

    ## samples with genotype data
    donors_cpg = sample_mapping['InternalID'].unique()
    donors_cpg.sort()
    donors_geno = sorted(set(list(G.sample.values)).intersection(donors_cpg))
    logging.info(
        'Number of unique donors with genotype data: {}'.format(len(donors_geno))
    )

    ## samples with both (can this be done in one step?)
    sample_mapping1 = sample_mapping.loc[sample_mapping['OneK1K_ID'].isin(donors_exprs)]
    sample_mapping_both = sample_mapping1.loc[
        sample_mapping1['InternalID'].isin(donors_geno)
    ]
    donors_e = sample_mapping_both['OneK1K_ID'].unique()
    donors_g = sample_mapping_both['InternalID'].unique()
    assert len(donors_e) == len(donors_g)

    logging.info('Number of unique common donors: {}'.format(len(donors_g)))

    ## samples in kinship
    donors_e_short = [re.sub('.*_', '', donor) for donor in donors_e]
    donors_k = sorted(set(list(K.sample_0.values)).intersection(donors_e_short))

    ######################################
    ### subset files to common samples ###
    ######################################

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
    G = G.sel(sample=donors_g)

    # make data frame to save as csv
    data = G.values
    Z_df = pd.DataFrame(data, columns=G.snp.values, index=G.sample.values)
    Z_df = Z_df.dropna(axis=1)

    # delete large files to free up memory
    del G

    ###############
    #### kinship
    K = K.sel(sample_0=donors_k, sample_1=donors_k)
    assert all(K.sample_0 == donors_k)
    assert all(K.sample_1 == donors_k)

    # make data frame to save as csv
    K_df = pd.DataFrame(K.values, columns=K.sample_0, index=K.sample_1)

    del K  # delete K to free up memory

    #########################################
    ############## save files ###############
    #########################################

    with expression_filename.open('w') as ef:
        y_df.to_csv(ef, index=False)

    with genotype_filename.open('w') as gf:
        Z_df.to_csv(gf, index=False)

    with kinship_filename.open('w') as kf:
        K_df.to_csv(kf, index=False)


if __name__ == '__main__':
    main()
