#!/bin/bash

function run() {
	local batch=$1
	local analysis_suf=$2
	local main_suf=$3
	local web_suf=$4
	local dir=work/${analysis_suf}

	test -f ${dir}/gender.tsv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/gender.tsv ${dir}/gender.tsv
	test -f ${dir}/age.csv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/age.csv ${dir}/age.csv
	test -f ${dir}/qc.csv || gsutil cp "gs://cpg-tob-wgs-${main_suf}/gvcf/batch${batch}/*.csv" ${dir}/qc.csv
	R --vanilla <<code
rmarkdown::render('qc.Rmd', output_file='qc.html', params=list(\
gender_tsv='${dir}/gender.tsv', \
age_csv='${dir}/age.csv', \
qc_csv='${dir}/qc.csv', \
gvcf_bucket_suffix='${main_suf}'\
))
code
	gsutil cp qc.html gs://cpg-tob-wgs-${web_suf}/qc/qc.html
}

# Run test first
run 0 test test temporary
