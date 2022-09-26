### Analysis plan

This folder contains three scripts:
* igll5_vep.py is a Python script which takes a MT + HT object, filters relevant variants (detail below) and exports to plink files
* prepare_inputs.py is a Python script which takes in genotype and expression files and prepares the input files to run SKAT
* run_SKAT.R is an R script that runs SKAT using the inputs generated from the previous script

#### Step 1 - select genetic variants

At the moment, I have subsetted both the MT object and the VEP-annotated object manually, using [this script]() to a specific genomic region around this one specific gene (IGLL5), in the future, either 1) add a previous step subsetting the MT + HT objects to the right genomic region taking a gene name as input, or add that to this script.

This step selects QC-passing, biallelic SNP that are rare (alternative allele frequency < 5%) and that are predicted by VEP to have regulatory consequences.
Then, it creates plink files (.bed, .bim, .fam) for those variants only.

#### Step 2 - prepare input files for SKAT

[SKAT]() is an R package to run gene-set association tests.
It requires as inputs:
* a genotype matrix (Z), samples X genetic variants
* a phenotype vector (y.c), samples X 1 (in my case, this will be the expression level of a gene across samples, see below)
* (optionally,) a kinship matrix (K), samples X samples
* (optionally,) a covariance matrix (X), samples X covariates

I use this (Python) script to prepare these input files, using similar steps as in the [CellRegMap pipeline]().

#### Step 3 - run SKAT

This is an R script that run SKAT-O in multiple modes:
* with and without the kinship matrix
* reporting the p-values for each of the burden, SKAT and optimised SKAT-O tests.