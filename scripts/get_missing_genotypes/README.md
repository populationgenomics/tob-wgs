# Investigate NA values in TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to get GT and GQ information for TOB samples. This is run on individuals with 'NA' values within the original TOB-WGS dataset and compared to GT calls and GQ scores within the DRAGEN version of the dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "tob_gt_gq/v0" \
--description "na GT" python3 main.py
```
