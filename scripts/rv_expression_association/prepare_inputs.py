#!/usr/bin/env python3

import logging
import pandas as pd
import scanpy as sc
import xarray as xr
from cloudpathlib import AnyPath
from pandas_plink import read_plink1_bin
from limix.qc import quantile_gaussianize

# use logging to print statements, display at info level
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)

phenotype_file = "gs://cpg-tob-wgs-test/v0/skat/sce22.h5ad"

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
    ######################################
    ###### sample mapping file (SMF) #####
    ######################################

    ## this file will map cells to donors
    # sample_mapping = pd.read_csv(
    #     sample_mapping_file,
    #     dtype={
    #         "individual_long": str,
    #         "genotype_individual_id": str,
    #         "phenotype_sample_id": str,
    #     },
    #     index_col=0,
    # )

    # ## extract unique individuals
    # donors0 = sample_mapping["genotype_individual_id"].unique()
    # donors0.sort()
    # logging.info("Number of unique donors: {}".format(len(donors0)))

    ######################################
    ########### phenotype file ###########
    ######################################

    #### TO DO: create pseudobulk

    # open anndata
    # parse the phenotype from the file (via write to temp)
    with open('i_am_a_temporary.h5ad', 'w', encoding='utf-8') as handle:
        handle.write(AnyPath(phenotype_file).read_text())
    adata = sc.read('i_am_a_temporary.h5ad')
    # adata = sc.read(phenotype_file)

    # sparse to dense
    mat = adata.raw.X.todense()
    # make pandas dataframe
    mat_df = pd.DataFrame(
        data=mat.T, index=adata.raw.var.index, columns=adata.obs.index
    )
    print(mar_df.head())
    # turn into xr array
    # phenotype = xr.DataArray(
    #     mat_df.values,
    #     dims=["trait", "cell"],
    #     coords={"trait": mat_df.index.values, "cell": mat_df.columns.values},
    # )
    # phenotype = phenotype.sel(cell=sample_mapping["phenotype_sample_id"].values)

    # # delete large files to free up memory
    # del mat
    # del mat_df

    # # select gene
    # y = phenotype.sel(trait=gene_name)
    # y = quantile_gaussianize(y)
    # y.c = y.values.reshape(y.shape[0], 1)

    # del phenotype  # delete to free up memory

    # ######################################
    # ############ kinship file ############
    # ######################################

    # ## read in GRM (genotype relationship matrix; kinship matrix)
    # K = pd.read_csv(kinship_file, index_col=0)
    # K.index = K.index.astype("str")
    # assert all(K.columns == K.index)  # symmetric matrix, donors x donors

    # K = xr.DataArray(
    #     K.values,
    #     dims=["sample_0", "sample_1"],
    #     coords={"sample_0": K.columns, "sample_1": K.index},
    # )
    # K = K.sortby("sample_0").sortby("sample_1")
    # donors = sorted(set(list(K.sample_0.values)).intersection(donors0))
    # logging.info("Number of donors after kinship intersection: {}".format(len(donors)))

    # ## subset to relevant donors
    # K = K.sel(sample_0=donors, sample_1=donors)
    # assert all(K.sample_0 == donors)
    # assert all(K.sample_1 == donors)

    # del K  # delete K to free up memory

    # ######################################
    # ############ genotype file ###########
    # ######################################

    # #### TO DO: create matrix from plink (no need to expand?)

    # ## read in genotype file (plink format)
    # G = read_plink1_bin(genotype_file)

    # ## select relavant SNPs based on feature variant filter file
    # fvf = pd.read_csv(feature_variant_file, index_col=0)
    # leads = fvf[fvf["feature"] == gene_name]["snp_id"].unique()
    # G_sel = G[:, G["snp"].isin(leads)]

    # # expand out genotypes from cells to donors (and select relevant donors in the same step)
    # G_expanded = G_sel.sel(sample=sample_mapping["individual_long"].values)
    # # assert all(hK_expanded.sample.values == G_expanded.sample.values)

    # # delete large files to free up memory
    # del G
    # del G_sel

    ### ouputs: y.c, Z, K (optional: X)

if __name__ == "__main__":
    main()