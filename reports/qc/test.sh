#!/bin/bash

function run() {
	local batch=$1
	local main_suf=$2
	local joint_calling_run_version=$3

	local dir=work/${main_suf}

	local test="FALSE"
	if [ $main_suf = "test" ]; then
		test="TRUE"
	fi

	test -f ${dir}/gender.tsv || gsutil cp "gs://cpg-tob-wgs-test/gender.tsv" ${dir}/gender.tsv
	test -f ${dir}/age.csv    || gsutil cp "gs://cpg-tob-wgs-test/age.csv" ${dir}/age.csv
	test -f ${dir}/qc.csv     || gsutil cp "gs://cpg-tob-wgs-${main_suf}/gvcf/${batch}/*.csv" ${dir}/qc.csv
	test -f ${dir}/meta.tsv   || gsutil cp "gs://cpg-tob-wgs-${main_suf}/joint-calling/${joint_calling_run_version}/meta.tsv"
	R --vanilla <<code
rmarkdown::render('qc.Rmd', output_file='qc.html', params=list(\
test=${test}, \
gender_tsv='${dir}/gender.tsv', \
age_csv='${dir}/age.csv', \
qc_csv='${dir}/qc.csv', \
meta_tsv='${dir}/meta.tsv', \
gvcf_bucket_suffix='${main_suf}'\
))
code
	gsutil cp qc.html gs://cpg-tob-wgs-${main_suf}-web/qc/qc-${joint_calling_run_version}.html
}

run batch1 test v2
