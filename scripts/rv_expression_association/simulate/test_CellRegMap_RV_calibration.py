#!/usr/bin/env python3
# pylint: disable=missing-module-docstring,wrong-import-position,import-error

# This script records differences in the calibration of gene-set associations
# under different scenarios. It is a mirrored script of test_CellRegMap.py,
# evaluating the same tests (CellRegMap-RV VT, burden and omnibus tests) and
# scenarios (1, 2, 2a, 3, 3a, 4), but critically shuffling the genotypes prior to
# testing, to assess the rate of false positives (calibration of tests)

# import python modules
import sys
import subprocess
import pandas as pd
import numpy as np
from cloudpathlib import AnyPath
from cpg_utils.hail_batch import output_path
from numpy import arange, eye, ones, zeros
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

from cellregmap import (
    run_gene_set_association,
    run_burden_association,
    omnibus_set_association,
)


def get_crm_pvs(pheno, covs, genotypes, E=None):
    """
    CellRegMap-RV tests
    * score test (variance)
    * burden test (max, sum, comphet)
    * omnibus (Cauchy) test
    """
    pv0 = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=E)[0]
    pv1 = run_burden_association(y=pheno, G=genotypes, W=covs, E=E, mask="mask.max")[0]
    pv2 = run_burden_association(y=pheno, G=genotypes, W=covs, E=E, mask="mask.sum")[0]
    pv3 = run_burden_association(
        y=pheno, G=genotypes, W=covs, E=E, mask="mask.comphet"
    )[0]
    pv4 = omnibus_set_association(np.array([pv0, pv1]))
    pv5 = omnibus_set_association(np.array([pv0, pv2]))
    pv6 = omnibus_set_association(np.array([pv0, pv3]))
    return [pv0, pv1, pv2, pv3, pv4, pv5, pv6]


# get genotypes
# these are variants in and around gene VPREB3 on chrom 22
g_file = AnyPath(output_path('VPREB3_50K_window/SNVs.csv'))
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
geno_subset = geno_all.sample(n_samples)
variant_count = geno_subset.sum(axis=0)  # get alt allele count
variant_freq = variant_count / (2 * n_samples)  # get alt allele frequency
# remove variants left all 0's after donor sub-sampling
variant_freq = variant_freq[variant_freq > 0]

# consider singletons (1 copy in 1 individual) only
singleton_freq = 0.5 / n_samples
all_singletons = list(variant_freq[variant_freq == singleton_freq].index.values)
print(len(all_singletons))

seed(0)
noise = randn(n_samples, 1)  # random noise Gaussian
seed(0)
noise_pois = poisson(lam=1, size=n_samples).reshape(n_samples, 1)  # Poisson noise
covs = ones((n_samples, 1))  # intercept of ones as covariates
E = eye(n_samples)

# same columns for all coming tables
cols = [
    'P_shapiro',
    'P_CRM_RV',
    'P_CRM_burden_max',
    'P_CRM_burden_sum',
    'P_CRM_burden_comphet',
    'P_CRM_omnibus_max',
    'P_CRM_omnibus_sum',
    'P_CRM_omnibus_comphet',
    'P_shapiro_Pois',
    'P_CRM_RV_Pois',
    'P_CRM_burden_max_Pois',
    'P_CRM_burden_sum_Pois',
    'P_CRM_burden_comphet_Pois',
    'P_CRM_omnibus_max_Pois',
    'P_CRM_omnibus_sum_Pois',
    'P_CRM_omnibus_comphet_Pois',
]

# scenario 1
# * test only those 10 variants
# * same direction and magnitude of effect
n_reps = 100
pv_scenario1_mt = zeros((n_reps, 16))
for i in range(n_reps):
    seed(i)
    singletons = all_singletons
    select_singletons_10 = sample(singletons, 10)
    genotypes = geno_subset[select_singletons_10]  # subset genotypes
    beta = ones((genotypes.shape[1], 1))  # create effect size
    # get other singletons to test (to assess calibration)
    [singletons.remove(x) for x in select_singletons_10]
    alt_singletons_10 = sample(list(singletons), 10)
    alt_genotypes = geno_subset[alt_singletons_10]  # subset genotypes
    # now build pheno with regular genotypes, test alt
    # Gaussian
    pheno = genotypes @ beta + noise  # build phenotype (Gauss)
    pv_scenario1_mt[i, 0] = shapiro(pheno).pvalue  # record normality pv
    pv_scenario1_mt[i, 1:8] = get_crm_pvs(pheno, covs, alt_genotypes, E)
    # Poisson
    pheno_pois = genotypes @ beta + noise_pois  # build phenotype Poisson
    pv_scenario1_mt[i, 8] = shapiro(pheno_pois).pvalue  # record normality pv
    pv_scenario1_mt[i, 9:17] = get_crm_pvs(pheno_pois, covs, alt_genotypes, E)

print(pv_scenario1_mt.shape)
print(len(cols))
print(len(['rep' + str(rep) for rep in range(n_reps)]))

pv_scenario1_df = pd.DataFrame(
    data=pv_scenario1_mt,
    columns=cols,
    index=['rep' + str(rep) for rep in range(n_reps)],
)

print(pv_scenario1_df.head())

pv_scenario1_filename = AnyPath(
    output_path('simulations/CRM/1000samples_10causal_singletons/shuffled/10tested_samebeta.csv')
)
with pv_scenario1_filename.open('w') as pf:
    pv_scenario1_df.to_csv(pf, index=False)


# scenario 2
# * test 50 variants (of which only 10 are causal)
# * same direction and magnitude of effects
pv_scenario2_mt = zeros((n_reps, 16))
for i in range(n_reps):
    seed(i)
    singletons = all_singletons
    select_singletons_50 = sample(singletons, 50)
    genotypes = geno_subset[select_singletons_50]  # subset genotypes
    beta = zeros((genotypes.shape[1], 1))        # create betas as 0s
    beta[0:10] = 1                               # only 10 non-0 betas
    # get other singletons to test (to assess calibration)
    [singletons.remove(x) for x in select_singletons_50]
    alt_singletons_50 = sample(singletons, 50)
    alt_genotypes = geno_subset[alt_singletons_50]  # subset genotypes
    # now build pheno with regular genotypes, test alt
    # Gaussian 
    pheno = genotypes @ beta + noise               # build phenotype (Gauss)
    pv_scenario2_mt[i, 0] = shapiro(pheno).pvalue  # record normality pv
    pv_scenario2_mt[i, 1:8] = get_crm_pvs(pheno, covs, alt_genotypes, E)
    # Poisson
    pheno_pois = genotypes @ beta + noise_pois          # build phenotype Poisson
    pv_scenario2_mt[i, 8] = shapiro(pheno_pois).pvalue  # record normality pv
    pv_scenario2_mt[i, 9:17] = get_crm_pvs(pheno_pois, covs, alt_genotypes, E)

pv_scenario2_df = pd.DataFrame(
    data=pv_scenario2_mt,
    columns=cols,
    index=['rep' + str(rep) for rep in range(n_reps)],
)

print(pv_scenario2_df.head())

pv_scenario2_filename = AnyPath(output_path('simulations/CRM/1000samples_10causal_singletons/shuffled/50tested_samebeta.csv'))
with pv_scenario2_filename.open('w') as pf:
    pv_scenario2_df.to_csv(pf, index=False)


# scenario 2a
# * test 20 variants (of which only 10 are causal)
# * same direction and magnitude of effects
pv_scenario2a_mt = zeros((n_reps, 16))
for i in range(n_reps):
    seed(i)
    singletons = all_singletons
    select_singletons_20 = sample(singletons, 20)
    genotypes = geno_subset[select_singletons_20]  # subset genotypes
    beta = zeros((genotypes.shape[1], 1))          # create betas as 0s
    beta[0:10] = 1                                 # only 10 non-0 betas
    # get other singletons to test (to assess calibration)
    [singletons.remove(x) for x in select_singletons_20]
    alt_singletons_20 = sample(singletons, 20)
    alt_genotypes = geno_subset[alt_singletons_20]  # subset genotypes
    # now build pheno with regular genotypes, test alt
    # Gaussian
    pheno = genotypes @ beta + noise                # build phenotype (Gauss)
    pv_scenario2a_mt[i, 0] = shapiro(pheno).pvalue  # record normality pv
    pv_scenario2a_mt[i, 1:8] = get_crm_pvs(pheno, covs, alt_genotypes, E)
    # Poisson
    pheno_pois = genotypes @ beta + noise_pois           # build phenotype Poisson
    pv_scenario2a_mt[i, 8] = shapiro(pheno_pois).pvalue  # record normality pv
    pv_scenario2a_mt[i, 9:17] = get_crm_pvs(pheno_pois, covs, alt_genotypes, E)

pv_scenario2a_df = pd.DataFrame(
    data=pv_scenario2a_mt,
    columns=cols,
    index=['rep' + str(rep) for rep in range(n_reps)],
)

print(pv_scenario2a_df.head())

pv_scenario2a_filename = AnyPath(output_path('simulations/CRM/1000samples_10causal_singletons/shuffled/20tested_samebeta.csv'))
with pv_scenario2a_filename.open('w') as pf:
    pv_scenario2a_df.to_csv(pf, index=False)


# scenario 3
# * test 10 variants
# * same magnitude of effect
# * vary direction for 2/10 variants
pv_scenario3_mt = zeros((n_reps, 16))
for i in range(n_reps):
    seed(i)
    singletons = all_singletons
    select_singletons_10 = sample(singletons, 10)
    genotypes = geno_subset[select_singletons_10]  # subset genotypes
    beta = ones((genotypes.shape[1], 1))           # create betas as 1s
    beta[0:2] = -1                                 # for two variants, -1
    # get other singletons to test (to assess calibration)
    [singletons.remove(x) for x in select_singletons_10]
    alt_singletons_10 = sample(singletons, 10)
    alt_genotypes = geno_subset[alt_singletons_10]  # subset genotypes
    # now build pheno with regular genotypes, test alt
    # Gaussian
    pheno = genotypes @ beta + noise               # build phenotype (Gauss)
    pv_scenario3_mt[i, 0] = shapiro(pheno).pvalue  # record normality pv
    pv_scenario3_mt[i, 1:8] = get_crm_pvs(pheno, covs, alt_genotypes, E)
    # Poisson
    pheno_pois = genotypes @ beta + noise_pois          # build phenotype Poisson
    pv_scenario3_mt[i, 8] = shapiro(pheno_pois).pvalue  # record normality pv
    pv_scenario3_mt[i, 9:17] = get_crm_pvs(pheno_pois, covs, alt_genotypes, E)

pv_scenario3_df = pd.DataFrame(
    data=pv_scenario3_mt,
    columns=cols,
    index=['rep' + str(rep) for rep in range(n_reps)],
)

print(pv_scenario3_df.head())

pv_scenario3_filename = AnyPath(output_path('simulations/CRM/1000samples_10causal_singletons/shuffled/10tested_2negativebeta.csv'))
with pv_scenario3_filename.open('w') as pf:
    pv_scenario3_df.to_csv(pf, index=False)


# scenario 3a
# * test 10 variants
# * same magnitude of effect
# * vary direction for 5/10 variants
pv_scenario3a_mt = zeros((n_reps, 16))
for i in range(n_reps):
    seed(i)
    singletons = all_singletons
    select_singletons_10 = sample(list(singletons), 10)
    genotypes = geno_subset[select_singletons_10]  # subset genotypes
    beta = ones((genotypes.shape[1], 1))           # create betas as 1s
    beta[0:5] = -1                                 # for five variants, -1
    # get other singletons to test (to assess calibration)
    [singletons.remove(x) for x in select_singletons_10]
    alt_singletons_10 = sample(singletons, 10)
    alt_genotypes = geno_subset[alt_singletons_10]  # subset genotypes
    # now build pheno with regular genotypes, test alt
    # Gaussian
    pheno = genotypes @ beta + noise                # build phenotype (Gauss)
    pv_scenario3a_mt[i, 0] = shapiro(pheno).pvalue  # record normality pv
    pv_scenario3a_mt[i, 1:8] = get_crm_pvs(pheno, covs, alt_genotypes, E)
    # Poisson
    pheno_pois = genotypes @ beta + noise_pois           # build phenotype Poisson
    pv_scenario3a_mt[i, 8] = shapiro(pheno_pois).pvalue  # record normality pv
    pv_scenario3a_mt[i, 9:17] = get_crm_pvs(pheno_pois, covs, alt_genotypes, E)

pv_scenario3a_df = pd.DataFrame(
    data=pv_scenario3a_mt,
    columns=cols,
    index=['rep' + str(rep) for rep in range(n_reps)],
)

print(pv_scenario3a_df.head())

pv_scenario3a_filename = AnyPath(output_path('simulations/CRM/1000samples_10causal_singletons/shuffled/10tested_5negativebeta.csv'))
with pv_scenario3a_filename.open('w') as pf:
    pv_scenario3a_df.to_csv(pf, index=False)


# scenario 4
# * test 10 variants
# * same direction of effect
# * vary magnitude
pv_scenario4_mt = zeros((n_reps, 16))
for i in range(n_reps):
    seed(i)
    singletons = all_singletons
    select_singletons_10 = sample(singletons, 10)
    genotypes = geno_subset[select_singletons_10]  # subset genotypes
    beta = arange(0.1, 1.1, 0.1).reshape(10, 1)    # create varying betas
    # get other singletons to test (to assess calibration)
    [singletons.remove(x) for x in select_singletons_10]
    alt_singletons_10 = sample(singletons, 10)
    alt_genotypes = geno_subset[alt_singletons_10]  # subset genotypes
    # now build pheno with regular genotypes, test alt
    # Gaussian
    pheno = genotypes @ beta + noise               # build phenotype (Gauss)
    pv_scenario4_mt[i, 0] = shapiro(pheno).pvalue  # record normality pv
    pv_scenario4_mt[i, 1:8] = get_crm_pvs(pheno, covs, alt_genotypes, E)
    # Poisson
    pheno_pois = genotypes @ beta + noise_pois          # build phenotype Poisson
    pv_scenario4_mt[i, 8] = shapiro(pheno_pois).pvalue  # record normality pv
    pv_scenario4_mt[i, 9:17] = get_crm_pvs(pheno_pois, covs, alt_genotypes, E)

pv_scenario4_df = pd.DataFrame(
    data=pv_scenario4_mt,
    columns=cols,
    index=['rep' + str(rep) for rep in range(n_reps)],
)

print(pv_scenario4_df.head())

pv_scenario4_filename = AnyPath(output_path('simulations/CRM/1000samples_10causal_singletons/shuffled/10tested_varyingbeta.csv'))
with pv_scenario4_filename.open('w') as pf:
    pv_scenario4_df.to_csv(pf, index=False)
