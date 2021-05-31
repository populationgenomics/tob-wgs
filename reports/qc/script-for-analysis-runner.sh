#!/bin/bash

set -ex

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -p|--prod)
    PROD=YES
    shift # past argument
    ;;
    -b|--batch)
    BATCH="$2"
    shift # past argument
    shift # past value
    ;;
    -v|--version)
    VERSION="$2"
    shift # past argument
    shift # past value
    ;;
esac
done

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
	r-reactable r-ggrepel r-sessioninfo r-gargle r-here r-assertthat r-dt r-googlecloudstorager r-ggforce

function run() {
	local batch=$1
	local analysis_suf=$2
	local main_suf=$3
	local joint_calling_run_version=$4
	local is_test=$5

	local dir=work/${analysis_suf}
	local test="FALSE"
	local meta_csv_path="gs://cpg-tob-wgs-temporary/joint-calling/${joint_calling_run_version}/sample_qc/meta.tsv"
	if [ $is_test -eq 1 ]; then
		test="TRUE"
		meta_csv_path="gs://cpg-tob-wgs-test/joint-calling/test-${joint_calling_run_version}/sample_qc/meta.tsv"
	fi

	test -f ${dir}/gender.tsv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/gender.tsv ${dir}/gender.tsv
	test -f ${dir}/age.csv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/age.csv ${dir}/age.csv
	test -f ${dir}/qc.csv || gsutil cp "gs://cpg-tob-wgs-${main_suf}/gvcf/batch${batch}/*.csv" ${dir}/qc.csv
	test -f ${dir}/meta.tsv || gsutil cp $meta_csv_path ${dir}/meta.tsv
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

# Run on full data in a standard access level
if [[ $PROD = "YES" ]]
then
	run $BATCH analysis main $VERSION 0
else
	run 1 test test $VERSION 1
fi
