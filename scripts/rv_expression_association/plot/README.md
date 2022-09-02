# Plot data using the analysis runner

<!-- This runs a Hail query script in Dataproc using Hail Batch in order to plot a test dataset using the analysis runner. -->
Remember to run ```gcloud auth application-default login``` first

```sh
analysis-runner --dataset tob-wgs --description "plot alternative allele frequencies" --output-dir "plot/v0" --access-level test python3 plot_alt_af.py
```

```sh
analysis-runner --dataset tob-wgs --description "plot qc metrics" --output-dir "plot/v0" --access-level test python3 plot_qc.py
```

```sh
analysis-runner --dataset tob-wgs --description "print number of variants at AF threshold" --output-dir "plot/v0" --access-level test python3 extract_rare_variants.py
```
