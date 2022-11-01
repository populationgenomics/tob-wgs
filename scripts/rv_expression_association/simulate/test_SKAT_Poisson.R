#!/usr/bin/env Rscript

## This script aims to record differences in the power of gene-set associations
## under different scenarios. It runs a SKAT, burden, and SKAT-O tests
## (implemented in the SKAT R package; Wu et al AJHG 2011, Lee et al AJHG 2012)
## to test for an association between a set of rare genetic variants (real)
## and a simulated phenotype. The sets of variants tested and the simulated
## effect sizes are varied and the resulted association p-values are recorded.

## specifically in this script I model the phenotype as Poisson distributed

# install R packages
install.packages("googleCloudStorageR", repos = "http://cran.csiro.au")
install.packages("SKAT", repos = "http://cran.csiro.au/")

print("Install STAAR dependencies Rcpp")

# required for STAAR
install.packages("Rcpp", repos = "http://cran.csiro.au")
install.packages("RcppArmadillo", repos = "http://cran.csiro.au")
library(Rcpp)
library(RcppArmadillo)

print("Install GENESIS dependencies igraph")
install.packages("igraph", repos = "http://cran.csiro.au")
library(igraph)

print("Install STAAR dependencies GENESIS")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GENESIS") # this currently fails because of igraph

print("Install STAAR dependencies SeqArray")
# required for GMMAT
BiocManager::install("SeqArray")
BiocManager::install("SeqVarTools")

print("Install STAAR dependencies GMMAT")
# install GMMAT
install.packages("GMMAT", repos = "http://cran.csiro.au")

print("Install STAAR")
# install STAAR
library(devtools)
devtools::install_github("xihaoli/STAAR")

# import R packages
library(googleCloudStorageR)
library(gargle)
library(SKAT)
library(STAAR)
library(glue)

pvalues <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
print(CCT(pvals = pvalues))

# # token authorisation (Google Cloud Storage R)
# scope <- c("https://www.googleapis.com/auth/cloud-platform")
# token <- token_fetch(scopes = scope)
# gcs_auth(token = token)

# # set bucket
# googleCloudStorageR::gcs_global_bucket("gs://cpg-tob-wgs-test")

# # get genotypes
# # these are variants in and around gene VPREB3 on chrom 22
# g_file <- googleCloudStorageR::gcs_get_object("v0/VPREB3_50K_window/SNVs.csv")
# g_df <- as.data.frame(g_file)

# # because the current matrix is counting the copies of the reference allele
# # while we are interested in the alternative allele, flip the genotypes
# geno_all <- 2 - as.matrix(g_df)

# n_samples <- 1000

# # Step 1: sample 1000 individuals only
# # first setting
# # * 1000 individuals
# # * 10 causal variants one in each if 10 individuals

# set.seed(0)
# geno_1000 <- geno_all[sample(nrow(geno_all), 1000), ]
# variant_count <- colSums(geno_1000)             # get alt allele count
# variant_freq <- variant_count / (2 * n_samples) # get alt allele frequency
# # remove variants left all 0"s after donor sub-sampling
# variant_freq <- variant_freq[variant_freq %in% variant_freq[variant_freq > 0]]

# # consider singletons (1 copy in 1 individual) only
# singleton_freq <- 0.5 / n_samples
# singletons <- names(variant_freq[variant_freq == singleton_freq])

# set.seed(1)
# rnoise <- rnorm(n_samples)
# print(shapiro.test(rnoise))
# covs <- matrix(1, nrow = n_samples, ncol = 1) # intercept of ones as covariates

# # scenario 1
# # * test only those 10 variants
# # * same direction and magnitude of effect
# n_reps <- 1000
# pv_scenario1_mt <- matrix(0, nrow = n_reps, ncol = 5)
# for (i in 1:n_reps){
#     set.seed(i)
#     select_singletons_10 <- singletons[sample(length(singletons), 10)]
#     genotypes <- geno_1000[, select_singletons_10]       # subset genotypes
#     beta <- matrix(1, nrow = ncol(genotypes), ncol = 1)  # create effect size
#     lambda <- genotypes %*% beta                         # get mean parameter
#     set.seed(0)
#     pheno <- rpois(n = n_samples, lambda = lambda)  # build phenotype (Poisson)
#     pheno <- pheno + rnoise                         # add random normal noise
#     pv_normal <- shapiro.test(pheno)$p.value             # record normality pv
#     obj <- SKAT_Null_Model(pheno ~ covs, out_type = "C") # build null model SKAT
#     pv_skat <- SKAT(genotypes, obj)$p.value                     # SKAT
#     pv_burden <- SKAT(genotypes, obj, r.corr = 1)$p.value       # burden
#     pv_skat_o <- SKAT(genotypes, obj, method = "SKATO")$p.value # SKAT-O
#     pv_scenario1_mt[i, 1] <- pv_normal
#     pv_scenario1_mt[i, 2] <- pv_skat
#     pv_scenario1_mt[i, 3] <- pv_burden
#     pv_scenario1_mt[i, 4] <- pv_skat_o
#     pv_scenario1_mt[i, 5] <- var(lambda) / (var(lambda) + 1)
# }
# pv_scenario1_df <- as.data.frame(pv_scenario1_mt)
# colnames(pv_scenario1_df) <- c("P_shapiro", "P_SKAT", "P_burden", "P_SKATO", "geno_beta_var")
# rownames(pv_scenario1_df) <- paste0("rep", 1:n_reps)

# print(head(pv_scenario1_df))

# pv_scenario1_filename <- "10tested_samebeta.csv"
# write.csv(pv_scenario1_df, pv_scenario1_filename)

# # scenario 2
# # * test 50 variants (of which only 10 are causal)
# # * same direction and magnitude of effect
# pv_scenario2_mt <- matrix(0, nrow = n_reps, ncol = 5)
# for (i in 1:n_reps){
#     set.seed(i)
#     select_singletons_50 <- singletons[sample(length(singletons), 50)]
#     genotypes <- geno_1000[, select_singletons_50]       # subset genotypes
#     beta <- matrix(0, nrow = ncol(genotypes), ncol = 1)  # create betas as 0s
#     beta[1:10] <- 1                                      # only 10 non-0 betas
#     lambda <- genotypes %*% beta                         # get mean parameter
#     set.seed(0)
#     pheno <- rpois(n = n_samples, lambda = lambda)  # build phenotype (Poisson)
#     pheno <- pheno + rnoise                         # add random normal noise
#     pv_normal <- shapiro.test(pheno)$p.value             # record normality pv
#     obj <- SKAT_Null_Model(pheno ~ covs, out_type = "C") # build null model SKAT
#     pv_skat <- SKAT(genotypes, obj)$p.value                     # SKAT
#     pv_burden <- SKAT(genotypes, obj, r.corr = 1)$p.value       # burden
#     pv_skat_o <- SKAT(genotypes, obj, method = "SKATO")$p.value # SKAT-O
#     pv_scenario2_mt[i, 1] <- pv_normal
#     pv_scenario2_mt[i, 2] <- pv_skat
#     pv_scenario2_mt[i, 3] <- pv_burden
#     pv_scenario2_mt[i, 4] <- pv_skat_o
#     pv_scenario2_mt[i, 5] <- var(lambda) / (var(lambda) + 1)
# }
# pv_scenario2_df <- as.data.frame(pv_scenario2_mt)
# colnames(pv_scenario2_df) <- c("P_shapiro", "P_SKAT", "P_burden", "P_SKATO", "geno_beta_var")
# rownames(pv_scenario2_df) <- paste0("rep", 1:n_reps)

# print(head(pv_scenario2_df))

# pv_scenario2_filename <- "50tested_samebeta.csv"
# write.csv(pv_scenario2_df, pv_scenario2_filename)

# # scenario 2a
# # * test 20 variants (of which only 10 are causal)
# # * same direction and magnitude of effect
# pv_scenario2a_mt <- matrix(0, nrow = n_reps, ncol = 5)
# for (i in 1:n_reps){
#     set.seed(i)
#     select_singletons_20 <- singletons[sample(length(singletons), 20)]
#     genotypes <- geno_1000[, select_singletons_20]       # subset genotypes
#     beta <- matrix(0, nrow = ncol(genotypes), ncol = 1)  # create betas as 0s
#     beta[1:10] <- 1                                      # only 10 non-0 betas
#     lambda <- genotypes %*% beta                         # get mean parameter
#     set.seed(0)
#     pheno <- rpois(n = n_samples, lambda = lambda)  # build phenotype (Poisson)
#     pheno <- pheno + rnoise                         # add random normal noise
#     pv_normal <- shapiro.test(pheno)$p.value             # record normality pv
#     obj <- SKAT_Null_Model(pheno ~ covs, out_type = "C") # build null model SKAT
#     pv_skat <- SKAT(genotypes, obj)$p.value                     # SKAT
#     pv_burden <- SKAT(genotypes, obj, r.corr = 1)$p.value       # burden
#     pv_skat_o <- SKAT(genotypes, obj, method = "SKATO")$p.value # SKAT-O
#     pv_scenario2a_mt[i, 1] <- pv_normal
#     pv_scenario2a_mt[i, 2] <- pv_skat
#     pv_scenario2a_mt[i, 3] <- pv_burden
#     pv_scenario2a_mt[i, 4] <- pv_skat_o
#     pv_scenario2a_mt[i, 5] <- var(lambda) / (var(lambda) + 1)
# }
# pv_scenario2a_df <- as.data.frame(pv_scenario2a_mt)
# colnames(pv_scenario2a_df) <- c("P_shapiro", "P_SKAT", "P_burden", "P_SKATO", "geno_beta_var")
# rownames(pv_scenario2a_df) <- paste0("rep", 1:n_reps)

# print(head(pv_scenario2a_df))

# pv_scenario2a_filename <- "20tested_samebeta.csv"
# write.csv(pv_scenario2a_df, pv_scenario2a_filename)

# # scenario 3
# # * test 10 variants
# # * same magnitude of effect
# # * vary direction for 2/10 variants
# pv_scenario3_mt <- matrix(0, nrow = n_reps, ncol = 5)
# for (i in 1:n_reps){
#     set.seed(i)
#     select_singletons_10 <- singletons[sample(length(singletons), 10)]
#     genotypes <- geno_1000[, select_singletons_10]       # subset genotypes
#     beta <- matrix(1, nrow = ncol(genotypes), ncol = 1)  # create betas as 1s
#     beta[1:2] <- -1                                      # for two variants, -1
#     lambda <- genotypes %*% beta                         # get mean parameter
#     set.seed(0)
#     pheno <- rpois(n = n_samples, lambda = lambda)  # build phenotype (Poisson)
#     pheno <- pheno + rnoise                         # add random normal noise
#     pv_normal <- shapiro.test(pheno)$p.value             # record normality pv
#     obj <- SKAT_Null_Model(pheno ~ covs, out_type = "C") # build null model SKAT
#     pv_skat <- SKAT(genotypes, obj)$p.value                     # SKAT
#     pv_burden <- SKAT(genotypes, obj, r.corr = 1)$p.value       # burden
#     pv_skat_o <- SKAT(genotypes, obj, method = "SKATO")$p.value # SKAT-O
#     pv_scenario3_mt[i, 1] <- pv_normal
#     pv_scenario3_mt[i, 2] <- pv_skat
#     pv_scenario3_mt[i, 3] <- pv_burden
#     pv_scenario3_mt[i, 4] <- pv_skat_o
#     pv_scenario3_mt[i, 5] <- var(lambda) / (var(lambda) + 1)
# }
# pv_scenario3_df <- as.data.frame(pv_scenario3_mt)
# colnames(pv_scenario3_df) <- c("P_shapiro", "P_SKAT", "P_burden", "P_SKATO", "geno_beta_var")
# rownames(pv_scenario3_df) <- paste0("rep", 1:n_reps)

# print(head(pv_scenario3_df))

# pv_scenario3_filename <- "10tested_2negativebeta.csv"
# write.csv(pv_scenario3_df, pv_scenario3_filename)

# # scenario 4
# # * test 10 variants
# # * same direction of effect
# # * vary magnitude
# pv_scenario4_mt <- matrix(0, nrow = n_reps, ncol = 5)
# for (i in 1:n_reps){
#     set.seed(i)
#     select_singletons_10 <- singletons[sample(length(singletons), 10)]
#     genotypes <- geno_1000[, select_singletons_10]       # subset genotypes
#     beta <- seq(0.1, 1, by = 0.1)                        # create varying betas
#     lambda <- genotypes %*% beta                         # get mean parameter
#     set.seed(0)
#     pheno <- rpois(n = n_samples, lambda = lambda)  # build phenotype (Poisson)
#     pheno <- pheno + rnoise                         # add random normal noise
#     pv_normal <- shapiro.test(pheno)$p.value             # record normality pv
#     obj <- SKAT_Null_Model(pheno ~ covs, out_type = "C") # build null model SKAT
#     pv_skat <- SKAT(genotypes, obj)$p.value                     # SKAT
#     pv_burden <- SKAT(genotypes, obj, r.corr = 1)$p.value       # burden
#     pv_skat_o <- SKAT(genotypes, obj, method = "SKATO")$p.value # SKAT-O
#     pv_scenario4_mt[i, 1] <- pv_normal
#     pv_scenario4_mt[i, 2] <- pv_skat
#     pv_scenario4_mt[i, 3] <- pv_burden
#     pv_scenario4_mt[i, 4] <- pv_skat_o
#     pv_scenario4_mt[i, 5] <- var(lambda) / (var(lambda) + 1)
# }
# pv_scenario4_df <- as.data.frame(pv_scenario4_mt)
# colnames(pv_scenario4_df) <- c("P_shapiro", "P_SKAT", "P_burden", "P_SKATO", "geno_beta_var")
# rownames(pv_scenario4_df) <- paste0("rep", 1:n_reps)

# print(head(pv_scenario4_df))

# pv_scenario4_filename <- "10tested_varyingbeta.csv"
# write.csv(pv_scenario4_df, pv_scenario4_filename)

# # attempt at saving using code from
# # https://github.com/populationgenomics/analysis-runner/blob/main/examples/r/script.R
# dataset_env <- Sys.getenv("tob-wgs")
# gcs_outdir <- glue("gs://cpg-tob-wgs-test/v0/simulations/skat/1000samples_10causal_singletons_poisson/")
# system(glue("gsutil cp {pv_scenario1_filename} {gcs_outdir}"))
# system(glue("gsutil cp {pv_scenario2_filename} {gcs_outdir}"))
# system(glue("gsutil cp {pv_scenario2a_filename} {gcs_outdir}"))
# system(glue("gsutil cp {pv_scenario3_filename} {gcs_outdir}"))
# system(glue("gsutil cp {pv_scenario4_filename} {gcs_outdir}"))
# cat(glue("[{date()}] Finished successfully!"))