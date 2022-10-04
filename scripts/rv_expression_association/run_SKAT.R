#!/usr/bin/env Rscript
""""
This script runs SKAT (Wu et al AJHG 2011, Lee et al AJHG 2012) to test for an association 
between a set of rare genetic variants and the expression level of a gene. 
First, it makes sure that the samples are matching and in the same order
across all objects.
""""

install.packages("googleCloudStorageR", repos = 'http://cran.csiro.au')
install.packages("SKAT", repos = 'http://cran.csiro.au/')

library(googleCloudStorageR)
library(gargle)
library(dplyr)
library(SKAT)

# token authorisation (Google Cloud Storage R)
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)

# set bucket
googleCloudStorageR::gcs_global_bucket("gs://cpg-tob-wgs-test")

# sample mapping file, matching CPG IDs to OneK1K IDs
sample_mapping_file <- googleCloudStorageR::gcs_get_object("scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv")
smf_df <- as.data.frame(sample_mapping_file)
colnames(smf_df)[2] <- "CPG_ID"  # change column name to match other files

# expression file (gene: VPREB3, cell type: naive B cells)
exprs <- googleCloudStorageR::gcs_get_object("v0/vpreb3_B_naive_expression.csv")
exprs_df <- as.data.frame(exprs)
colnames(exprs_df)[1] <- 'OneK1K_ID'  # change column name to match other files

# genotype file (rare variants, freq<5%, predicted by VEP to have regulatory consequences)
geno <- googleCloudStorageR::gcs_get_object("v0/vpreb3_rare_regulatory.csv")
geno_df <- as.data.frame(geno)
colnames(geno_df)[1] <- "CPG_ID"  # again match column names

# expression files use the OneK1K IDs, genotypes the CPG IDs, use inner join and the SMF to match
tmp_df <- inner_join(smf_df, exprs_df)
print(paste0("number of samples with expression data: ", nrow(tmp_df)))
tmp_df1 <- inner_join(tmp_df, geno_df)
print(paste0("number of samples with expression and genotype data: ", nrow(tmp_df1)))

# extract phenotype values (.c for continuous, using SKAT notation)
y.c <- matrix(tmp_df1[,4], ncol = 1)
print(paste0("y.c dimensionality: ", dim(y.c)))

# extract genotype values (Z, using SKAT notation)
Z <- as.matrix(tmp_df1[, 5:ncol(tmp_df1)])
print(paste0("Z dimensionality: ", dim(Z)))

# genotypes are inverted (figure out!)
Z <- 2-Z

# add covariate matrix (X, using SKAT notation)
X <- matrix(1, nrow = nrow(tmp_df1), ncol = 1)   # intercept only
print(paste0("X dimensionality: ", dim(X)))

# run null model (no Kinship)
obj <- SKAT_Null_Model(y.c ~ X, out_type = "C")

# get all three p-values
pv_skat <- SKAT(Z, obj)$p.value                    # SKAT
pv_burden <- SKAT(Z, obj, r.corr = 1)$p.value      # burden
pv_skat_o <- SKAT(Z, obj, method="SKATO")$p.value  # SKAT-O

print(paste0("no Kinship p-values, SKAT: ", pv_skat, ", burden: ", pv_burden, ", SKAT-O: ", pv_skat_o))

###################################################
######### Adding Kinship matrix to the model

# kinship matrix uses yet another ID (a shorter version of the OneK1K ID)
# use this file I had generated previously to match the two
smf2 <- googleCloudStorageR::gcs_get_object("v0/skat/smf_Bcells.csv")
smf2_df <- as.data.frame(smf2)
colnames(smf2_df)[3:4] <- c("OneK1K_short_ID","OneK1K_ID")
smf2_df$OneK1K_short_ID <- as.character(smf2_df$OneK1K_short_ID)  # modify type
smf2_df$OneK1K_ID <- as.character(smf2_df$OneK1K_ID)
smf2_df <- smf2_df[, c("OneK1K_short_ID", "OneK1K_ID")]  # do not consider cell IDs for these purposes
smf2_df <- smf2_df[!duplicated(smf2_df), ]  # remove duplicate rows

# merge sample mapping files
tmp_df2 <- inner_join(tmp_df1, smf2_df)

# load kinship file
Kin <- googleCloudStorageR::gcs_get_object("v0/skat/grm_wide.csv")
K_df <- as.data.frame(Kin)
colnames(K_df)[1] <- "OneK1K_short_ID"

# select relevant samples
tmp_df3 <- inner_join(tmp_df2[, c("OneK1K_ID","OneK1K_short_ID")], K_df)
K <- as.matrix(tmp_df3[, 3:ncol(tmp_df3)])
K <- K[, tmp_df3$OneK1K_short_ID]
print(paste0("K dimensionality: ", dim(K)))

# run null model (with Kinship matrix)
obj_K <- SKAT_NULL_emmaX(y.c ~ X, K = K)

# get all three p-values
pv_skat <- SKAT(Z, obj_K)$p.value  # SKAT
pv_burden <- SKAT(Z, obj_K, r.corr = 1)$p.value  # burden
pv_skat_o <- SKAT(Z, obj_K, method = "SKATO")$p.value  # SKAT-O

print(paste0("with Kinship p-values, SKAT: ", pv_skat, ", burden: ", pv_burden, ", SKAT-O: ", pv_skat_o))
