#!/usr/bin/env python3

## This script aims to record differences in the power of gene-set associations
## under different scenarios. It runs our CellRegMap-RV method
## to test for an association between a set of rare genetic variants (real)
## and a simulated phenotype. The sets of variants tested and the simulated
## effect sizes are varied and the resulted association p-values are recorded.

# import python modules
import sys
import subprocess
import pandas as pd
from cloudpathlib import AnyPath
from cpg_utils.hail_batch import output_path
from numpy import eye, ones, zeros
from random import sample
from numpy.random import poisson, randn, seed
from scipy.stats import shapiro

# install CellRegMap (new version) from github
subprocess.run(
    [
        sys.executable,
        '-m',
        'pip',
        'install',
        'git+https://github.com/annacuomo/CellRegMap',
        '--force-reinstall',  # install github version and overwrite current
        '--no-dependencies',  # same dependencies, no need to uninstall and reinstall those
    ],
    check=True,
)

from cellregmap import run_gene_set_association


# get genotypes
# these are variants in and around gene VPREB3 on chrom 22
g_file = AnyPath(output_path("VPREB3_50K_window/SNVs.csv"))
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
variant_count = geno_1000.sum(axis=0)  # get alt allele count
variant_freq = variant_count / (2 * n_samples)  # get alt allele frequency
# remove variants left all 0"s after donor sub-sampling
variant_freq = variant_freq[variant_freq > 0]

# consider singletons (1 copy in 1 individual) only
singleton_freq = 0.5 / n_samples
singletons = variant_freq[variant_freq == singleton_freq].index.values

seed(0)
noise = randn(n_samples, 1)  # random noise Gaussian
noise_pois = poisson(lam=1, size=n_samples).reshape(n_samples, 1)  # Poisson noise
covs = ones((n_samples, 1))  # intercept of ones as covariates
E = eye(n_samples)

# scenario 1
# * test only those 10 variants
# * same direction and magnitude of effect
n_reps = 1000
pv_scenario1_mt = zeros((n_reps, 3))
for i in range(n_reps):
    seed(i)
    select_singletons_10 = sample(list(singletons), 10)
    genotypes = geno_1000[select_singletons_10]  # subset genotypes
    beta = ones((genotypes.shape[1], 1))  # create effect size
    pheno = genotypes @ beta + noise  # build phenotype
    pheno_pois = genotypes @ beta + noise_pois  # build phenotype Poisson
    pv_normal = shapiro(pheno).pvalue  # record normality pv
    pv_crm_rv = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=E)[0]
    pv_crm_rv_pois = run_gene_set_association(y=pheno_pois, G=genotypes, W=covs, E=E)[0]
    pv_scenario1_mt[i, 0] = pv_normal
    pv_scenario1_mt[i, 1] = pv_crm_rv
    pv_scenario1_mt[i, 2] = pv_crm_rv_pois

pv_scenario1_df = pd.DataFrame(
    data=pv_scenario1_mt,
    columns=["P_shapiro", "P_CRM_RV", "P_CRM_RV_Pois"],
    index=["rep" + str(rep) for rep in range(n_reps)],
)

print(pv_scenario1_df.head())

pv_scenario1_filename = AnyPath(output_path("simulations/CRM/1000samples_10causal_singletons/10tested_samebeta.csv"))
with pv_scenario1_filename.open('w') as pf:
    pv_scenario1_df.to_csv(pf, index=False)


# scenario 2
# * test 50 variants (of which only 10 are causal)
# * same direction and magnitude of effects
pv_scenario2_mt = zeros((n_reps, 3))
for i in range(n_reps):
    seed(i)
    select_singletons_50 = sample(list(singletons), 50)
    genotypes = geno_1000[select_singletons_50] # subset genotypes
    beta = zeros((genotypes.shape[1], 1))       # create betas as 0s
    betas[0:10] = 1                             # only 10 non-0 betas
    pheno = genotypes @ beta + noise            # build phenotype Gauss
    pheno_pois = genotypes @ beta + noise_pois  # build phenotype Poisson
    pv_normal = shapiro(pheno).pvalue  # record normality pv
    pv_crm_rv = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=E)[0]
    pv_crm_rv_pois = run_gene_set_association(y=pheno_pois, G=genotypes, W=covs, E=E)[0]
    pv_scenario2_mt[i, 0] = pv_normal
    pv_scenario2_mt[i, 1] = pv_crm_rv
    pv_scenario2_mt[i, 2] = pv_crm_rv_pois

pv_scenario2_df = pd.DataFrame(
    data=pv_scenario2_mt,
    columns=["P_shapiro", "P_CRM_RV", "P_CRM_RV_Pois"],
    index=["rep" + str(rep) for rep in range(n_reps)],
)

print(pv_scenario2_df.head())

pv_scenario2_filename = AnyPath(output_path("simulations/CRM/1000samples_10causal_singletons/50tested_samebeta.csv"))
with pv_scenario2_filename.open('w') as pf:
    pv_scenario2_df.to_csv(pf, index=False)


# scenario 2a
# * test 20 variants (of which only 10 are causal)
# * same direction and magnitude of effects
pv_scenario2a_mt = zeros((n_reps, 3))
for i in range(n_reps):
    seed(i)
    select_singletons_20 = sample(list(singletons), 20)
    genotypes = geno_1000[select_singletons_20] # subset genotypes
    beta = zeros((genotypes.shape[1], 1))       # create betas as 0s
    betas[0:10] = 1                             # only 10 non-0 betas
    pheno = genotypes @ beta + noise            # build phenotype Gauss
    pheno_pois = genotypes @ beta + noise_pois  # build phenotype Poisson
    pv_normal = shapiro(pheno).pvalue           # record normality pv
    pv_crm_rv = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=E)[0]
    pv_crm_rv_pois = run_gene_set_association(y=pheno_pois, G=genotypes, W=covs, E=E)[0]
    pv_scenario2a_mt[i, 0] = pv_normal
    pv_scenario2a_mt[i, 1] = pv_crm_rv
    pv_scenario2a_mt[i, 2] = pv_crm_rv_pois

pv_scenario2a_df = pd.DataFrame(
    data=pv_scenario2a_mt,
    columns=["P_shapiro", "P_CRM_RV", "P_CRM_RV_Pois"],
    index=["rep" + str(rep) for rep in range(n_reps)],
)

print(pv_scenario2a_df.head())

pv_scenario2a_filename = AnyPath(output_path("simulations/CRM/1000samples_10causal_singletons/20tested_samebeta.csv"))
with pv_scenario2a_filename.open('w') as pf:
    pv_scenario2a_df.to_csv(pf, index=False)

