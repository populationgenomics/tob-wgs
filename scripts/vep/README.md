# Annotate with VEP using the analysis-runner

This runs a Hail query script in Dataproc using Hail Batch in order to run VEP on a Hail matrix table. To run, install the analysis-runner, then execute the following command:

```sh
VEP_VERSION=108.2 \
    analysis-runner \
    --dataset tob-wgs \
    --description "Annotate with VEP $VEP_VERSION" \
    --output-dir tob_wgs_vep/v7_vep_$VEP_VERSION \
    --access-level standard \
    main.py \
    --script run_vep.py \
    --mt gs://cpg-tob-wgs-main/mt/v7.mt \
    --vep-version $VEP_VERSION
```
