#!/bin/bash

function run() {
	local batch=$1
	local analysis_suf=$2
	local main_suf=$3
	local web_suf=$4
	local joint_calling_run_version=$5
	local is_test=$6

	local dir=work/${analysis_suf}
	local test="FALSE"
	local out_version="${joint_calling_run_version}"
	local meta_csv_path="gs://cpg-tob-wgs-temporary/joint-calling/${joint_calling_run_version}/sample_qc/meta.tsv"
	if [ $is_test -eq 1 ]; then
		test="TRUE"
		out_version="test-${out_version}"
		meta_csv_path="gs://cpg-tob-wgs-test/joint-calling/test-${joint_calling_run_version}/sample_qc/meta.tsv"
	fi

	test -f ${dir}/gender.tsv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/gender.tsv ${dir}/gender.tsv
	test -f ${dir}/age.csv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/age.csv ${dir}/age.csv
	test -f ${dir}/qc.csv || gsutil cp "gs://cpg-tob-wgs-${main_suf}/gvcf/batch${batch}/*.csv" ${dir}/qc.csv
	test -f ${dir}/meta.tsv || gsutil cp gs://cpg-tob-wgs-temporary/joint-calling/${joint_calling_run_version}/sample_qc/meta.tsv ${dir}/meta.tsv
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
	gsutil cp qc.html gs://cpg-tob-wgs-web/qc/qc-${out_version}.html
}


# Run test first
run 1 test test temporary v1 1
