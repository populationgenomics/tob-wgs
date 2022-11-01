#!/usr/bin/env Rscript

## This script aims to record differences in the power of gene-set associations
## under different scenarios. It runs a SKAT, burden, and SKAT-O tests
## (implemented in the SKAT R package; Wu et al AJHG 2011, Lee et al AJHG 2012)
## to test for an association between a set of rare genetic variants (real)
## and a simulated phenotype. The sets of variants tested and the simulated
## effect sizes are varied and the resulted association p-values are recorded.

## specifically in this script I model the phenotype as Poisson distributed

# install R packages (install.packages)
install.packages("googleCloudStorageR", repos = "http://cran.csiro.au")
install.packages("SKAT", repos = "http://cran.csiro.au/")
# required for STAAR
install.packages("Rcpp", repos = "http://cran.csiro.au")
install.packages("RcppArmadillo", repos = "http://cran.csiro.au")
# Bioconductor packages install
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.csiro.au")
BiocManager::install("GENESIS")     # required for STAAR
BiocManager::install("SeqArray")    # required for GMMAT
BiocManager::install("SeqVarTools") # required for GMMAT
install.packages("GMMAT", repos = "http://cran.csiro.au") # STAAR dep.

# install STAAR and ACAT using devtools
library(devtools)
devtools::install_github("yaowuliu/ACAT")
devtools::install_github("xihaoli/STAAR")

# import R packages
library(googleCloudStorageR)
library(gargle)
library(SKAT)
library(STAAR)
library(glue)

# token authorisation (Google Cloud Storage R)
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)

# set bucket
googleCloudStorageR::gcs_global_bucket("gs://cpg-tob-wgs-test")

# SKAT tests (SKAT, burden, SKAT-O)
# Wu et al AJHG 2011, Lee et al AJHG 2012
get_skat_pvs <- function(pheno, covs, genotypes, weights = c(1, 25)) {
    obj <- SKAT::SKAT_Null_Model(pheno ~ covs, out_type = "C")  # null
    pv_skat <- SKAT::SKAT(genotypes, obj,                       # SKAT
        weights.beta = weights)$p.value 
    pv_burden <- SKAT::SKAT(genotypes, obj, r.corr = 1,         # burden
        weights.beta = weights)$p.value
    pv_skat_o <- SKAT::SKAT(genotypes, obj, method = "SKATO",   # SKAT-O
        weights.beta = weights)$p.value
    return(c(pv_skat, pv_burden, pv_skat_o))
}

# ACAT-V test (Liu et al, AJHG 2019)
get_acatv_pv <- function(pheno, covs, genotypes, weights = c(1, 25)) {
    obj_acat <- ACAT::NULL_Model(t(pheno), covs)                    # null
    pv <- ACAT::ACAT_V(genotypes, obj_acat, weights.beta = weights) # ACAT-V
    return(pv)
}

# ACAT-O omnibus test (Liu et al, AJHG 2019)
# takes p-values from various tests (SKAT, burden and ACAT-V) and
# different sets of params (1,1) and (1,25) and returns a combined p-value
get_acato_pv <- function(pvals) {
    elems <- tan((0.5 - pvals) * pi)
    t_acato <- (1 / length(pvals)) * sum(elems)            # T statistic
    pv <- as.numeric(pcauchy(t_acato, lower.tail = FALSE)) # get Cauchy PV
    return(pv)
}

# STAAR combined test (Li et al, Nature Methods 2022)
get_cct_pv <- function(pvals) {
    pv <- STAAR::CCT(pvals = pvals)
    return(pv)
}

# STAAR actual test (STAAR-O)
get_staar_pv <- function(pheno, covs, genotypes) {
    fixed <- t(pheno) ~ covs
    data <- data.frame(pheno = t(pheno), covs = covs)
    print(head(data))
    obj_null_model <- STAAR::fit_null_glm(fixed, data, family = gaussian)
    res <- STAAR::STAAR(genotypes, obj_null_model)
    pv <- res["results_STAAR_O"][1]
    return(pv)
}

# utility function to get p-values from the different tests
# from pheno, covs and genos
get_all_pvs <- function(pheno, covs, genotypes, n_tests) {
    # pvals <- as.vector(matrix(0, nrow = n_tests))
    # pvals[1] <- shapiro.test(pheno)$p.value          # record normality pv
    # # SKAT
    # pvals[2] <- get_skat_pvs(pheno, covs, genotypes, weights = c(1, 1))[1]
    # pvals[3] <- get_skat_pvs(pheno, covs, genotypes, weights = c(1, 25))[1]
    # # burden
    # pvals[4] <- get_skat_pvs(pheno, covs, genotypes, weights = c(1, 1))[2]
    # pvals[5] <- get_skat_pvs(pheno, covs, genotypes, weights = c(1, 25))[2]
    # # ACAT-V
    # pvals[6] <- get_acatv_pv(pheno, covs, genotypes, weights = c(1, 1))
    # pvals[7] <- get_acatv_pv(pheno, covs, genotypes, weights = c(1, 25))
    # # SKAT-O
    # pvals[8] <- get_skat_pvs(pheno, covs, genotypes, weights = c(1, 1))[3]
    # pvals[9] <- get_skat_pvs(pheno, covs, genotypes, weights = c(1, 25))[3]
    # # ACAT-O (combining SKAT, burden and ACAT-V)
    # pvals[10] <- get_acato_pv(pvals[2:7])
    # STAAR (combined CCT)
    # pvals[11] <- get_cct_pv(pvals[2:7])
    # STAAR-O
    pvals[12] <- get_staar_pv(pheno, covs, genotypes)
    return(pvals)
}

# get genotypes
# these are variants in and around gene VPREB3 on chrom 22
g_file <- googleCloudStorageR::gcs_get_object("v0/VPREB3_50K_window/SNVs.csv")
g_df <- as.data.frame(g_file)

# because the current matrix is counting the copies of the reference allele
# while we are interested in the alternative allele, flip the genotypes
geno_all <- 2 - as.matrix(g_df)

n_samples <- 100

# Step 1: sample 1000 individuals
# first setting
# * 1000 individuals
# * 10 causal variants one in each if 10 individuals

set.seed(0)
geno_1000 <- geno_all[sample(nrow(geno_all), 1000), ]
variant_count <- colSums(geno_1000)             # get alt allele count
variant_freq <- variant_count / (2 * n_samples) # get alt allele frequency
# remove variants left all 0"s after donor sub-sampling
variant_freq <- variant_freq[variant_freq %in% variant_freq[variant_freq > 0]]

# consider singletons (1 copy in 1 individual) only
singleton_freq <- 0.5 / n_samples
singletons <- names(variant_freq[variant_freq == singleton_freq])

set.seed(0)
noise <- rnorm(n_samples)                       # random noise (Gauss)
noise_pois <- rpois(n = n_samples, lambda = 1)  # random noise (Poisson)
covs <- matrix(rnorm(n_samples * 2), ncol = 2)  # random covariates

cols <- c("P_shapiro", "P_SKAT_1_1", "P_SKAT_1_25",
    "P_burden_1_1", "P_burden_1_25", "P_ACATV_1_1", "P_ACATV_1_25",
    "P_SKATO_1_1", "P_SKATO_1_25", "P_ACATO", "P_STAAR", "P_STAARO")
    # "P_shapiro_Pois", "P_SKAT_1_1_Pois", "P_SKAT_1_25_Pois",
    # "P_burden_1_1_Pois", "P_burden_1_25_Pois", "P_ACATV_1_1_Pois",
    # "P_ACATV_1_25_Pois", "P_SKATO_1_1_Pois", "P_SKATO_1_25_Pois",
    # "P_ACATO_Pois", "P_STAAR_Pois")

# scenario 1
# * test only those 10 variants
# * same direction and magnitude of effect
n_reps <- 10
pv_scenario1_mt <- matrix(0, nrow = n_reps, ncol = 12)
for (i in 1:n_reps){
    set.seed(i)
    select_singletons_10 <- singletons[sample(length(singletons), 10)]
    genotypes <- geno_1000[, select_singletons_10]       # subset genotypes
    beta <- matrix(1, nrow = ncol(genotypes), ncol = 1)  # create effect size
    # Gaussian noise
    pheno <- genotypes %*% beta + noise               # build phenotype (Gauss)
    pv_scenario1_mt[i, 1:12] <- get_all_pvs(pheno, covs, genotypes, 12)
    # Poisson noise
    # pheno_pois <- genotypes %*% beta + noise_pois     # build phenotype (Pois)
    # pv_scenario1_mt[i, 12:22] <- get_all_pvs(pheno_pois, covs, genotypes, 11)
}
pv_scenario1_df <- as.data.frame(pv_scenario1_mt)
colnames(pv_scenario1_df) <- cols
rownames(pv_scenario1_df) <- paste0("rep", 1:n_reps)

print(head(pv_scenario1_df))

pv_scenario1_filename <- "10tested_samebeta.csv"
write.csv(pv_scenario1_df, pv_scenario1_filename)


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