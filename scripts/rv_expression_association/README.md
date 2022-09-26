## Analysis plan

This folder contains three scripts:
* [igll5_vep.py](igll5_vep.py) is a Python script which takes a MT + HT object, filters relevant variants (detail below) and exports to plink files
* [prepare_inputs.py](prepare_inputs.py) is a Python script which takes in genotype and expression files and prepares the input files to run SKAT
* [run_SKAT.R](run_SKAT.R) is an R script that runs SKAT using the inputs generated from the previous script

### Step 1 - select genetic variants

At the moment, I have subsetted both the MT object and the VEP-annotated object manually, using [this script](https://github.com/populationgenomics/analysis-runner/blob/main/scripts/subset_matrix_table.py) to a specific genomic region around this one specific gene (IGLL5), e.g., for the VEP-annotated hail table:
```
analysis-runner --dataset tob-wgs --description "subset vep annotated ht" --output-dir "v0" --access-level standard python3 subset_hail_table.py -i gs://cpg-tob-wgs-main/tob_wgs_vep/104/vep104.3_GRCh38.ht --chr chr22 --pos 22837780-22946111 --out IGLL5_50K_window_vep
```
in the future, either 1) add a previous step subsetting the MT + HT objects to the right genomic region taking a gene name as input, or add that to this script.

This step selects QC-passing, biallelic SNP that are rare (alternative allele frequency < 5%) and that are predicted by VEP to have regulatory consequences.
Then, it creates plink files (.bed, .bim, .fam) for those variants only.
To run this, use:
```
```

### Step 2 - prepare input files for SKAT

[SKAT]() is an R package to run gene-set association tests.
It requires as inputs:
* a genotype matrix (Z), samples X genetic variants
* a phenotype vector (y.c), samples X 1 (in my case, this will be the expression level of a gene across samples, see below)
* (optionally,) a kinship matrix (K), samples X samples
* (optionally,) a covariance matrix (X), samples X covariates

I use this (Python) script to prepare these input files, using similar steps as in the [CellRegMap pipeline](https://github.com/populationgenomics/cellregmap-pipeline).
To run:
```
analysis-runner --dataset "tob-wgs" \
    --description "open adata using scanpy" \
    --access-level "test" \
    --output-dir "tob_wgs_rv/expression_outliers" \
    prepare_inputs.py
```

### Step 3 - run SKAT

This is an R script that run SKAT-O in multiple modes:
* with and without the kinship matrix
* reporting the p-values for each of the burden, SKAT and optimised SKAT-O tests.
