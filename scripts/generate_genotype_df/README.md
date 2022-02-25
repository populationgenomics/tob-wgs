# Export TOB-WGS dataset to PLINK format

This runs a Hail query script in Dataproc using Hail Batch in order to export the TOB-WGS genotypes, which are requirements for the association analysis. Only variants that pass VQSR and GQ filters are retained, as well as variants with an MAF > 0.01. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_plink_maf01/v0" \
--description "TOB plink" python3 main.py
```
