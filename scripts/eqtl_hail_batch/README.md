# Run hail batch job to generate list of eQTLs

This runs a Hail batch script in order to generate a list of eQTLs from scRNA-seq expression. This [code](https://github.com/powellgenomicslab/onek1k_phase1/tree/main/single_cell_cis_eQTL_mapping) was taken from Seyhan Yazar from Joseph Powell's group at the Garvan-Weizmann Centre for Cellular Genomics, then converted into Python/hail batch. In contrast to the above script which was filtered for variants with a MAF <= 0.05, this pipeline filters the TOB-WGS dataset on the following criteria:
* Variants with MAF <= 0.01
* Variants which don’t pass the VQSR filter
* Low quality genotypes, set as those <= 20
* Samples with a genotype call rate < 0.8

## Workflow:

**Pre analysis**

Input: Cohort joint-call, map of sample identifiers, frequently table path, vep annotation path
Output: A filtered matrix table with rows (alleles) filtered on the following requirements:

1. biallelic
2. meets VQSR filters
3. gene quality score higher than 20,
4. call rate of 0.8, and 
5. variants with MAF <= 0.01.

**Base analysis (Round 1):**

Input: all SNP genotypes (run on each chromosome, full scRNA-seq expression data (run for each cell type)

Output: Set of significantly associated SNPs, expression data adjusted for covariates

Steps:
1. Fit a linear model to the expression data in order to account for experimental covariates (sex, age, PCs 1-4, Peer Factors 1 & 2). Take the residuals, which are then used as the adjusted expression values.
2. Test the association of all SNPs to the adjusted expression values

**Conditional analysis (Rounds 2-5):**

Input: significant SNP genotypes from the previous round, expression data adjusted for experimental covariates (or in round 3+, adjusted for experimental covariates and previous round lead SNPs)

Output: set of significantly associated SNPs, expression data adjusted for lead SNP & covariates

Steps:
1. Fit a model to the adjusted expression data to account for the genotype of the most significant SNP. Take the residuals, which are the SNP + covariate adjusted expression values
2. Test the association of significant SNPs from the previous round against SNP & covariate adjusted expression data

## Running via the analysis-runner:

To run, use conda to install the analysis-runner. For the first round of eQTLs, execute the following command for running one cell type and chromosome (B_intermediate cells and chromosome 22 are shown here as an example):

> Note, there is no way to run individual steps anymore, the single pipeline runs both eqtl + conditional analysis and uses checkpoints within the output-directory to avoid re-running stages

You can omit the `--chromosomes` and `--cell-types` parameters to run on all known chromosomes and cell_types from within the input-files-prefix directory

```sh
analysis-runner \
    --dataset tob-wgs \
    --access-level test \
    --output-dir "scrna-seq/eqtl_output/v3" \
    --description "eqtl batch job" \
    python3 launch_eqtl_spearman.py \
      --input-files-prefix gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files \
      --chromosomes 22 \
      --cell-types B_intermediate
```
