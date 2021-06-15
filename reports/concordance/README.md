# Concordance of SNPchip with WGS genotype calls

Here we're exploring genotype concordance between SNPchip and WGS data for the
TOB-WGS project. We have been provided with one multi-sample VCF file containing
genotype calls from the SNPchip data, and several single-sample GVCF files
containing genotype calls from the WGS data (merged into a Hail
[MatrixTable](https://hail.is/docs/0.2/hail.MatrixTable.html) (MT) for
processing).

## Contents

This directory consists of:

- [concordance.Rmd](concordance.Rmd): script written in R and Python which runs
  Hail Query on the SNPchip and WGS data, and generates a TSV file with the raw
  results and a HTML file displaying the summarised results.
- [concordance](concordance): CLI interface for running the above.
- [Dockerfile](Dockerfile): script to build the Docker image containing the two
  scripts above, running inside a conda environment containing Hail and several
  R packages.
- [environment.yml](environment.yml): contains conda packages available inside
  the Docker container.
- [batch_concordance.py](batch_concordance.py): Python script which submits a
  concordance job to Hail Batch.

In its current form, the Hail Query pipeline runs on a single VM using the
number of CPUs specified in the Batch script. Using a 16 CPU machine for testing
a MT containing 19 WGS samples had a wall clock time of 1 hour (for a chr22
subset, this went down to 3 minutes). This setup is okay for testing. For
running on larger datasets, Hail's
[dataproc](https://hail.is/docs/0.2/cloud/google_cloud.html#hailctl-dataproc)
setup would be a solution.

## Analysis Runner

In order to run this pipeline using the
[analysis runner](https://github.com/populationgenomics/analysis-runner),
you can use the following BASH command (make sure you've activated a conda
environment containing the
[analysis-runner](https://anaconda.org/cpg/analysis-runner) conda package):

```bash
analysis-runner \
  --access-level "test" \
  --dataset "tob-wgs" \
  --description "concordance" \
  --output-dir "gs://cpg-tob-wgs-test/snpchip/v1" \
  batch_concordance.py
```
