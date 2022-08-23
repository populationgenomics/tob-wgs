# Plot data using the analysis runner

This runs a Hail query script in Dataproc using Hail Batch in order to plot a test dataset using the analysis runner.
Remember to run ```sh gcloud auth application-default login``` first

```sh
analysis-runner --dataset tob-wgs --description "plot data" --output-dir "plot/v0" --access-level test python3 plot_main.py
```
