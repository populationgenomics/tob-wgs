# Run hail batch job to generate list of eQTLs

This runs a Hail batch script in order to generate a list of eQTLs from scRNA-seq expression. This code was taken from Seyhan Yazar from Joseph Powell's group at the Garvan-Weizmann Centre for Cellular Genomics, then converted into Python/hail batch. To run, use conda to install the analysis-runner. For the first round of eQTLs, execute the following command for running one cell type and chromosome:

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "kat/v0" \
    --description "eqtl batch job" \
    python3 generate_eqtl_spearman.py \
        --output-prefix 'gs://cpg-tob-wgs-test/kat/plasma/chr22/'
        --expression 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/expression_files/B_intermediate_expression.tsv' \
        --genotype 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/genotype_files/tob_genotype_chr22.tsv' \
        --geneloc 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv' \
        --snploc 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/snp_location_files/snpsloc_chr22.tsv' \
        --covariates 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/B_intermediate_peer_factors.tsv' \
        --keys 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "kat/v0" \
    --description "eqtl batch job" \
    python3 round2.conditional_analysis_test.py \
        --output-prefix 'gs://cpg-tob-wgs-test/kat/plasma/chr22/' \
        --residuals 'gs://cpg-tob-wgs-test/kat/plasma/chr22/log_residuals.tsv' \
        --significant-snps 'gs://cpg-tob-wgs-test/kat/plasma/chr22/correlation_results.csv' \
        --genotype 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/genotype_files/tob_genotype_chr22.tsv' \
        --keys 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv' \
        --test_subset_genes 5 # test with 5 genes only
```

To launch all cell types and chromosomes at once, run the following python wrapper script for the first round of eQTL analysis:

```sh
python3 launch_generate_eqtl_spearman.py \
--input-path "gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files" \
--output-dir 'gs://cpg-tob-wgs-test/eqtl_output' --chromosomes '22'
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
python3 launch_round2.conditional_analysis_test.py \
--input-path "gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files" \
--output-dir 'gs://cpg-tob-wgs-test/eqtl_output' --chromosomes '22' --first-round-path 'gs://cpg-tob-wgs-test/kat/'
```
