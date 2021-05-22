#!/bin/bash

test -f qc.csv || gstuil cp gs://cpg-tob-wgs-main/gvcf/batch2/R_210315_BINKAN1_1K1KDNA_M003.csv qc.csv
test -f gender.tsv || gsutil cp gs://cpg-tob-wgs-analysis/gender.tsv gender.tsv
test -f age.csv || gsutil cp gs://cpg-tob-wgs-analysis/age.csv age.csv

Rscript -e "rmarkdown::render('qc_exploration.Rmd', output_file='qc_exploration.html', params=list(gender_tsv='gender.tsv', age_csv='age.csv', qc_csv='qc.csv', bucket_suffix='test'))"

gsutil cp qc_exploration.html $OUTPUT/qc_exploration.html
