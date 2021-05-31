# Perform sample QC on the combined HGDP + 1KG and TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_sample_qc/v0" \
--description "hgdp1kg tobwgs sample qc" python3 main.py
```
