# Filter TOB-WGS dataset and export to parquet format

This runs a Hail query script in Dataproc using Hail Batch in order to export the TOB-WGS genotypes, which are requirements for the association analysis. Only variants that pass VQSR and GQ filters are retained, as well as variants with an MAF > 0.01. To run, use pip to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "scrna-seq/grch38_association_files/genotype_files/" \
--description "TOB parquet" python3 generate_genotype_df.py
```
