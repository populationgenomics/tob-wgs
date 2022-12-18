# Annotate with VEP using the analysis-runner

This runs a Hail query script in Dataproc using Hail Batch in order to run VEP on a Hail matrix table. TOB-WGS uses GENCODE 42, which corresponds to VEP 108. To generate annotated tables for both the `test` and `main` namespaces, install the analysis-runner, then execute the following command:

```sh
for ACCESS_LEVEL in test standard; do
    VEP_VERSION=108.2 \
    analysis-runner \
    --dataset tob-wgs \
    --description "Annotate with VEP $VEP_VERSION" \
    --output-dir tob_wgs_vep/v7_vep_$VEP_VERSION \
    --access-level $ACCESS_LEVEL \
    main.py \
    --script run_vep.py \
    --mt mt/v7.mt \
    --vep-version $VEP_VERSION
done
```
