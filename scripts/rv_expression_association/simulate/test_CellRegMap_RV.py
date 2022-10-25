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
from numpy import eye, ones, zeros
from random import sample
from numpy.random import poisson, randn, seed
from scipy.stats import shapiro

# install CellRegMap (new version) from github
subprocess.run([sys.executable,'-m','pip','install','git+https://github.com/annacuomo/CellRegMap'], check=True)

from cellregmap import run_gene_set_association


# get genotypes
# these are variants in and around gene VPREB3 on chrom 22
g_file = AnyPath(output_path("v0/VPREB3_50K_window/SNVs.csv"))
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
geno_1000 = geno_all.sample(n_samples)
variant_count = geno_1000.sum(axis=0)             # get alt allele count
variant_freq = variant_count / (2 * n_samples) # get alt allele frequency
# remove variants left all 0"s after donor sub-sampling
variant_freq = variant_freq[variant_freq > 0]

# consider singletons (1 copy in 1 individual) only
singleton_freq = 0.5 / n_samples
singletons = variant_freq[variant_freq == singleton_freq].index.values

seed(0)
noise = randn(n_samples, 1)            # random noise Gaussian
noise_pois = poisson(lam=1, size=n_samples).reshape(n_samples,1) # Poisson noise
covs = ones((n_samples, 1)) # intercept of ones as covariates
E = eye(n_samples)

# scenario 1
# * test only those 10 variants
# * same direction and magnitude of effect
n_reps = 1000
pv_scenario1_mt = zeros((n_reps, 2))
for i in range(n_reps):
    seed(i)
    select_singletons_10 = sample(list(singletons), 10)
    genotypes = geno_1000[select_singletons_10]       # subset genotypes
    beta = ones((genotypes.shape[1], 1))        # create effect size
    pheno = genotypes @ beta + noise            # build phenotype
    pheno_pois = genotypes @ beta + noise_pois  # build phenotype Poisson
    pv_normal = shapiro(pheno).pvalue           # record normality pv
    pv_crm_rv = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=E)[0] 
    pv_crm_rv_pois = run_gene_set_association(y=pheno_pois, G=genotypes, W=covs, E=E)[0]
    pv_scenario1_mt[i, 0] = pv_normal
    pv_scenario1_mt[i, 1] = pv_crm_rv
    pv_scenario1_mt[i, 2] = pv_crm_rv_pois

pv_scenario1_df = pd.DataFrame(data = pv_scenario1_mt, columns = ["P_shapiro", "P_CRM_RV", "P_CRM_RV_Pois"], index = "rep"+range(n_reps))
# )
# pv_scenario1_df.columns = ["P_shapiro", "P_CRM_RV", "P_CRM_RV_Pois")
# rownames(pv_scenario1_df) = "rep" + range (n_reps)

print(pv_scenario1_df.head())

pv_scenario1_filename = "10tested_samebeta.csv"
pv_scenario1_df.to_csv(pv_scenario1_filename)