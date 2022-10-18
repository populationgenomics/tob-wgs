#!/usr/bin/env Rscript

## This script aims to record differences in the power of gene-set associations
## under different scenarios. It runs a SKAT, burden, and SKAT-O tests
## (implemented in the SKAT R package; Wu et al AJHG 2011, Lee et al AJHG 2012)
## to test for an association between a set of rare genetic variants (real)
## and a simulated phenotype. The sets of variants tested and the simulated
## effect sizes are varied and the resulted association p-values are recorded.

# install R packages
install.packages("googleCloudStorageR", repos = 'http://cran.csiro.au')
install.packages("SKAT", repos = 'http://cran.csiro.au/')

# import R packages
library(googleCloudStorageR)
library(gargle)
library(SKAT)
library(glue)

# token authorisation (Google Cloud Storage R)
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)

# set bucket
googleCloudStorageR::gcs_global_bucket("gs://cpg-tob-wgs-test")
# set seed
set.seed(0)

# get genotypes
# these are variants in and around gene VPREB3 on chrom 22
G_file <- googleCloudStorageR::gcs_get_object("v0/VPREB3_50K_window/SNVs.csv")
G_df = as.data.frame(G_file)

# because the current matrix is counting the copies of the reference allele
# while we are interested in the alternative allele, flip the genotypes
Z = 2-as.matrix(G_df)

# Step 1: sample 100 individuals only
Z_100 = Z[sample(nrow(Z), 100), ]
variant_count = colSums(Z_100)       # get alt allele count
variant_freq = variant_count / 200   # get alt allele frequency
# remove variants left all 0's after donor sub-sampling
variant_freq = variant_freq[variant_freq %in% variant_freq[variant_freq>0]]

# consider singletons (1 copy in 1 individual) only
singletons = names(variant_freq[variant_freq==0.005])

# first scenario
# * 100 individuals
# * 10 variants one in each if 10 individuals
# * test only those 10 variants
# * same direction and magnitude of effect

n_samples = 100
noise = rnorm(n_samples)              # random noise
X = matrix(1, nrow=n_samples, ncol=1) # intercept of ones as covariates

n_reps = 100
pv_scenario1_mt = matrix(0, nrow=n_reps, ncol=4)
for (i in 1:n_reps){
    set.seed(i)
    select_singletons_10 = singletons[sample(length(singletons), 10)]
    G <- Z_100[,select_singletons_10]                   # subset genotypes
    beta = matrix(1, nrow=ncol(G), ncol=1)              # create effect size
    y = G %*% beta + noise                              # build phenotype
    pv_normal = shapiro.test(y)$p.value                 # record normality pv
    obj <- SKAT_Null_Model(y ~ X, out_type = "C")       # build null model SKAT
    pv_skat <- SKAT(G, obj)$p.value                     # SKAT
    pv_burden <- SKAT(G, obj, r.corr = 1)$p.value       # burden
    pv_skat_o <- SKAT(G, obj, method="SKATO")$p.value   # SKAT-O
    pv_scenario1_mt[i, 1] = pv_normal
    pv_scenario1_mt[i, 2] = pv_skat
    pv_scenario1_mt[i, 3] = pv_burden
    pv_scenario1_mt[i, 4] = pv_skat_o
}
pv_scenario1_df = as.data.frame(pv_scenario1_mt)
colnames(pv_scenario1_df) = c("P_shapiro","P_SKAT","P_burden","P_SKATO")
rownames(pv_scenario1_df) = paste0("rep",1:n_reps)

print(head(pv_scenario1_df))

pv_scenario1_filename = "10tested_samebeta.csv"
write.csv(pv_scenario1_df, pv_scenario1_filename)

# attempt at saving using code from
# https://github.com/populationgenomics/analysis-runner/blob/main/examples/r/script.R
dataset_env <- Sys.getenv("tob-wgs")
output_env <- Sys.getenv("v0/simulations/skat/100samples_10causalvariants/")
gcs_outdir <- glue("gs://cpg-tob-wgs-test/{output_env}")
system(glue("gsutil cp {pv_scenario1_filename} {gcs_outdir}"))
cat(glue("[{date()}] Finished successfully!"))