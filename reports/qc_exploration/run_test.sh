#!/bin/bash

test -f qc.csv || gsutil cp gs://cpg-tob-wgs-test/gvcf/batch0/R_210315_BINKAN1_1K1KDNA_M001.csv qc.csv
test -f gender.tsv || gsutil cp gs://cpg-tob-wgs-test/gender.tsv gender.tsv
test -f age.csv || gsutil cp gs://cpg-tob-wgs-test/age.csv age.csv

mamba install -y -c conda-forge r-reactable r-ggrepel r-sessioninfo r-gargle r-here r-assertthat r-dt r-googlecloudstorager
Rscript -e "rmarkdown::render('qc_exploration.Rmd', output_file='qc_exploration.html', params=list(gender_tsv='gender.tsv', age_csv='age.csv', qc_csv='qc.csv', bucket_suffix='test'))"

gsutil cp qc_exploration.html $OUTPUT/qc_exploration.html
