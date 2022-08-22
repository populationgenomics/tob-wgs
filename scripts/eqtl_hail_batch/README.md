# Run hail batch job to generate list of eQTLs

This runs a Hail batch script in order to generate a list of eQTLs from scRNA-seq expression. This [code](https://github.com/powellgenomicslab/onek1k_phase1/tree/main/single_cell_cis_eQTL_mapping) was taken from Seyhan Yazar from Joseph Powell's group at the Garvan-Weizmann Centre for Cellular Genomics, then converted into Python/hail batch. 

## Workflow:

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

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "scrna-seq/plasma/chr22/v0" \
    --description "eqtl batch job" \
    python3 generate_eqtl_spearman.py \
        --expression 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/expression_files/B_intermediate_expression.tsv' \
        --geneloc 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv' \
        --covariates 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/B_intermediate_peer_factors_file.txt' \
        --keys 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "scrna-seq/plasma/chr22/v0" \
    --description "eqtl batch job" \
    python3 conditional_analysis.py \
        --output-prefix 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v0' \
        --residuals 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v0/log_residuals.csv' \
        --significant-snps 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v0/correlation_results.tsv' \
        --test-subset-genes 5 # test with 5 genes only
```

To launch all cell types and chromosomes at once, run the following python wrapper script for the first round of eQTL analysis:

```sh
analysis-runner --dataset tob-wgs --access-level test --output-dir 'scrna-seq/eqtl_output/v0' --description "eqtl batch job" python3 launch_generate_eqtl_spearman.py --input-path "gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files" --chromosomes '22'
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
analysis-runner --dataset tob-wgs --access-level test --output-dir 'scrna-seq/eqtl_output/v0' --description "eqtl batch job" python3 launch_conditional_analysis.py --input-path "gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files" --chromosomes '22' --first-round-path 'gs://cpg-tob-wgs-test/scrna-seq/eqtl_output/v0' --output-dir 'gs://cpg-tob-wgs-test/scrna-seq/eqtl_output/v0'
```
