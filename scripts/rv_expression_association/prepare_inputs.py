#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import xarray as xr
from pandas_plink import read_plink1_bin

######################################
########### phenotype file ###########
######################################

# open anndata
adata = sc.read(phenotype_file)
# sparse to dense
mat = adata.raw.X.todense()
# make pandas dataframe
mat_df = pd.DataFrame(
    data=mat.T, index=adata.raw.var.index, columns=adata.obs.index
)
# turn into xr array
phenotype = xr.DataArray(
    mat_df.values,
    dims=["trait", "cell"],
    coords={"trait": mat_df.index.values, "cell": mat_df.columns.values},
)
phenotype = phenotype.sel(cell=sample_mapping["phenotype_sample_id"].values)

# delete large files to free up memory
del mat
del mat_df

# select gene
y = phenotype.sel(trait=gene_name)

y = quantile_gaussianize(y)
y.c = y.values.reshape(y.shape[0], 1)

del phenotype  # delete to free up memory

######################################
############ kinship file ############
######################################

## read in GRM (genotype relationship matrix; kinship matrix)
K = pd.read_csv(kinship_file, index_col=0)
K.index = K.index.astype("str")
assert all(K.columns == K.index)  # symmetric matrix, donors x donors

K = xr.DataArray(
    K.values,
    dims=["sample_0", "sample_1"],
    coords={"sample_0": K.columns, "sample_1": K.index},
)
K = K.sortby("sample_0").sortby("sample_1")
donors = sorted(set(list(K.sample_0.values)).intersection(donors0))
logging.info("Number of donors after kinship intersection: {}".format(len(donors)))

## subset to relevant donors
K = K.sel(sample_0=donors, sample_1=donors)
assert all(K.sample_0 == donors)
assert all(K.sample_1 == donors)

del K  # delete K to free up memory

######################################
############ genotype file ###########
######################################

## read in genotype file (plink format)
G = read_plink1_bin(genotype_file)

## select relavant SNPs based on feature variant filter file
fvf = pd.read_csv(feature_variant_file, index_col=0)
leads = fvf[fvf["feature"] == gene_name]["snp_id"].unique()
G_sel = G[:, G["snp"].isin(leads)]

# expand out genotypes from cells to donors (and select relevant donors in the same step)
G_expanded = G_sel.sel(sample=sample_mapping["individual_long"].values)
# assert all(hK_expanded.sample.values == G_expanded.sample.values)

# delete large files to free up memory
del G
del G_sel