#!/usr/bin/env python3

## This script aims to record differences in the power of gene-set associations
## under different scenarios. It runs our CellRegMap-RV method
## to test for an association between a set of rare genetic variants (real)
## and a simulated phenotype. The sets of variants tested and the simulated
## effect sizes are varied and the resulted association p-values are recorded.

# import python modules
import sys
import subprocesss
import pandas as pd
from cloudpathlib import AnyPath
from cpg_utils.hail_batch import output_path
from numpy import ones
from numpy.random import randn, sample, seed
from scipy.stats import shapiro

# install CellRegMap (new version) from github
subprocess.run([sys.executable,'-m','pip','install','git+https://github.com/annacuomo/CellRegMap'], check=True)

from cellregmap import run_gene_set_association


# get genotypes
# these are variants in and around gene VPREB3 on chrom 22
g_file = AnyPath(output_path(("v0/VPREB3_50K_window/SNVs.csv"))
g_df = pd.read_csv(g_file)

# because the current matrix is counting the copies of the reference allele
# while we are interested in the alternative allele, flip the genotypes
geno_all = 2 - g_df

n_samples = 1000

# Step 1: sample 1000 individuals only
# first setting
# * 1000 individuals
# * 10 causal variants one in each if 10 individuals

seed(0)
geno_1000 = geno_all[sample(nrow(geno_all), 1000), ]
variant_count = colSums(geno_1000)             # get alt allele count
variant_freq = variant_count / (2 * n_samples) # get alt allele frequency
# remove variants left all 0"s after donor sub-sampling
variant_freq = variant_freq[variant_freq %in% variant_freq[variant_freq > 0]]

# consider singletons (1 copy in 1 individual) only
singleton_freq = 0.5 / n_samples
singletons = names(variant_freq[variant_freq == singleton_freq])

seed(0)
noise = randn(n_samples, 1)                    # random noise
covs = ones((n_samples, 1)) # intercept of ones as covariates

# scenario 1
# * test only those 10 variants
# * same direction and magnitude of effect
n_reps = 1000
pv_scenario1_mt = matrix(0, nrow = n_reps, ncol = 2)
for i in range(n_reps):
    seed(i)
    select_singletons_10 = singletons[sample(length(singletons), 10)]
    genotypes = geno_1000[, select_singletons_10]       # subset genotypes
    beta = matrix(1, nrow = ncol(genotypes), ncol = 1)  # create effect size
    pheno = genotypes @ beta + noise                  # build phenotype
    pv_normal = shapiro(pheno).pvalue             # record normality pv
    pv_crm_rv = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=None) # TO DO allow E=None
    # obj = SKAT_Null_Model(pheno ~ covs, out_type = "C") # build null model SKAT
    # pv_skat = SKAT(genotypes, obj)$p.value                     # SKAT
    # pv_burden = SKAT(genotypes, obj, r.corr = 1)$p.value       # burden
    # pv_skat_o = SKAT(genotypes, obj, method = "SKATO")$p.value # SKAT-O
    pv_scenario1_mt[i, 1] = pv_normal
    pv_scenario1_mt[i, 2] = pv_crm_rv

pv_scenario1_df = pd.DataFrame(pv_scenario1_mt)
colnames(pv_scenario1_df) = c("P_shapiro", "P_CRM_RV")
rownames(pv_scenario1_df) = "rep" + range (n_reps)

print(head(pv_scenario1_df))

pv_scenario1_filename = "10tested_samebeta.csv"
write.csv(pv_scenario1_df, pv_scenario1_filename)