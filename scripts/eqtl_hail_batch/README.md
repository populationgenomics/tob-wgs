# Run hail batch job to generate list of eQTLs

This runs a Hail batch script in order to generate a list of eQTLs from scRNA-seq expression. This [code](https://github.com/powellgenomicslab/onek1k_phase1/tree/main/single_cell_cis_eQTL_mapping) was taken from Seyhan Yazar from Joseph Powell's group at the Garvan-Weizmann Centre for Cellular Genomics, then converted into Python/hail batch. To run, use conda to install the analysis-runner. For the first round of eQTLs, execute the following command for running one cell type and chromosome (B_intermediate cells and chromosome 22 are shown here as an example):

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "plasma/chr22/v0" \
    --description "eqtl batch job" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16 \
    python3 generate_eqtl_spearman.py \
        --output-prefix 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v0'
        --expression 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/expression_files/B_intermediate_expression.tsv' \
        --geneloc 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv' \
        --covariates 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/B_intermediate_peer_factors_file.txt' \
        --keys 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "plasma/chr22/v0" --output-dir "plasma/chr22/v0" \
    --description "eqtl batch job" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16 \
    python3 conditional_analysis.py \
        --output-prefix 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v5' \
        --residuals 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v5residual_df.tsv' \
        --significant-snps 'gs://cpg-tob-wgs-test/scrna-seq/plasma/chr22/v5correlation_results.csv' \
        --keys 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv' \
        --test-subset-genes 5 # test with 5 genes only
```

To launch all cell types and chromosomes at once, run the following python wrapper script for the first round of eQTL analysis:

```sh
python3 launch_generate_eqtl_spearman.py \
--input-path "gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files" \
--output-dir 'gs://cpg-tob-wgs-test/eqtl_output' --chromosomes '22'
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
python3 launch_conditional_analysis.py \
--input-path "gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files" \
--output-dir 'gs://cpg-tob-wgs-test/eqtl_output' --chromosomes '22' --first-round-path 'gs://cpg-tob-wgs-test/eqtl_output/'
```
