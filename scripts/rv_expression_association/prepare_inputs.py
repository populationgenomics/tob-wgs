#!/usr/bin/env python3

import subprocess
import sys

subprocess.run([sys.executable, '-m', 'pip', 'install', 'scanpy==1.9.1'], check=True)
subprocess.run([sys.executable, '-m', 'pip', 'install', 'limix'], check=True)
subprocess.run([sys.executable, '-m', 'pip', 'install', 'pandas_plink==2.2.9'], check=True)

import logging
import pandas as pd
import scanpy as sc
import xarray as xr
# from cloudpathlib import AnyPath
from cpg_utils import to_path
from pandas_plink import read_plink1_bin
from limix.qc import quantile_gaussianize

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# sample_mapping_file = 'gs://cpg-tob-wgs-test/v0/skat/smf_Bcells.csv'
phenotype_file = 'gs://cpg-tob-wgs-test/v0/skat/sce22.h5ad'
# genotype_file_bed = 'gs://cpg-tob-wgs-test/v0/skat/plink_chr22.bed'
# genotype_file_bim = 'gs://cpg-tob-wgs-test/v0/skat/plink_chr22.bim'
# genotype_file_fam = 'gs://cpg-tob-wgs-test/v0/skat/plink_chr22.fam'
# kinship_file = 'gs://cpg-tob-wgs-test/v0/skat/grm_wide.csv'


def main(
    # chrom: str,
    # gene_name: str,
    # sample_mapping_file: str,
    # genotype_file: str,
    # phenotype_file: str,
    # # context_file: str,
    # kinship_file: str,
    # # feature_variant_file: str,
    # output_folder: str,
    # # n_contexts: int = 10,
):
    #####################################
    ##### sample mapping file (SMF) #####
    #####################################

    # # this file will map cells to donors
    # sample_mapping = pd.read_csv(
    #     sample_mapping_file,
    #     dtype={
    #         'individual_long': str,
    #         'genotype_individual_id': str,
    #         'phenotype_sample_id': str,
    #     },
    #     index_col=0,
    # )

    # ## extract unique individuals
    # donors0 = sample_mapping['genotype_individual_id'].unique()
    # donors0.sort()
    # logging.info('Number of unique donors: {}'.format(len(donors0)))

    ######################################
    ########### phenotype file ###########
    ######################################

    #### TO DO: create pseudobulk

    # open anndata
    with to_path(phenotype_file).open('rb') as handle:
        data = handle.readlines()
    with open('temp.h5ad', 'wb') as handle:
        handle.writelines(data)
    adata = sc.read('temp.h5ad')

    # sparse to dense
    mat = adata.raw.X.todense()

    # create pseudobulk
    pseudobulk = mat[cells[0]].X.mean(axis=0)

    # make pandas dataframe
    mat_df = pd.DataFrame(
        data=mat.T, index=adata.raw.var.index, columns=adata.obs.index
    )
    print(mat_df.head())
    # turn into xr array
    phenotype = xr.DataArray(
        mat_df.values,
        dims=['trait', 'cell'],
        coords={'trait': mat_df.index.values, 'cell': mat_df.columns.values},
    )
    # phenotype = phenotype.sel(cell=sample_mapping['phenotype_sample_id'].values)

    # delete large files to free up memory
    del mat
    del mat_df

    # select gene
    y = phenotype.sel(trait=gene_name)
    y = quantile_gaussianize(y)
    y.c = y.values.reshape(y.shape[0], 1)

    del phenotype  # delete to free up memory

    # ######################################
    # ############ kinship file ############
    # ######################################

    # ## read in GRM (genotype relationship matrix; kinship matrix)
    # K = pd.read_csv(kinship_file, index_col=0)
    # K.index = K.index.astype('str')
    # assert all(K.columns == K.index)  # symmetric matrix, donors x donors

    # K = xr.DataArray(
    #     K.values,
    #     dims=['sample_0', 'sample_1'],
    #     coords={'sample_0': K.columns, 'sample_1': K.index},
    # )
    # K = K.sortby('sample_0').sortby('sample_1')
    # donors = sorted(set(list(K.sample_0.values)).intersection(donors0))
    # logging.info('Number of donors after kinship intersection: {}'.format(len(donors)))

    # ## subset to relevant donors
    # K = K.sel(sample_0=donors, sample_1=donors)
    # assert all(K.sample_0 == donors)
    # assert all(K.sample_1 == donors)

    # del K  # delete K to free up memory

    # ######################################
    # ############ genotype file ###########
    # ######################################

    #### TO DO: create matrix from plink (no need to expand?)

    ## read in genotype file (plink format)
    # bed
    # with to_path(genotype_file_bed).open('rb') as handle:
    #     data = handle.readlines()
    # with open('temp.bed', 'wb') as handle:
    #     handle.writelines(data)
    # # bim
    # with to_path(genotype_file_bim).open('rb') as handle:
    #     data = handle.readlines()
    # with open('temp.bim', 'wb') as handle:
    #     handle.writelines(data)
    # # fam
    # with to_path(genotype_file_fam).open('rb') as handle:
    #     data = handle.readlines()
    # with open('temp.fam', 'wb') as handle:
    #     handle.writelines(data)
    # # read
    # G = read_plink1_bin('temp.bed')
    
    # ## select relavant SNPs based on feature variant filter file
    # fvf = pd.read_csv(feature_variant_file, index_col=0)
    # set_variants = fvf[fvf['feature'] == gene_name]['snp_id'].unique()
    # G_sel = G[:, G['snp'].isin(set_variants)]

    # Z = G_sel.sel(sample=sample_mapping['individual_long'].values)
    # # expand out genotypes from cells to donors (and select relevant donors in the same step)
    # # G_expanded = G_sel.sel(sample=sample_mapping['individual_long'].values)
    # # assert all(hK_expanded.sample.values == G_expanded.sample.values)

    # # delete large files to free up memory
    # del G
    # del G_sel

    ### ouputs: y.c, Z, K (optional: X)

if __name__ == '__main__':
    main()