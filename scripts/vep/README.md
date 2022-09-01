# Test VEP using the analysis runner

This runs a Hail query script in Dataproc using Hail Batch in order to run VEP on a hail matrix table. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner \
	--dataset tob-wgs \
	--description "run vep" \
	--output-dir "tob_wgs_vep/v0" \
	--access-level main main.py \
	--script run_vep.py \
	--mt 'gs://cpg-tob-wgs-main/mt/v7.mt/' \
	--vep-version '104.3'
```
