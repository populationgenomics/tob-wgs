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
esac
done

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
	r-reactable r-ggrepel r-sessioninfo r-gargle r-here r-assertthat r-dt r-googlecloudstorager

function run() {
	local batch=$1
	local analysis_suf=$2
	local main_suf=$3
	local web_suf=$4

	test -f work/gender.tsv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/gender.tsv work/gender.tsv
	test -f work/age.csv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/age.csv work/age.csv
	test -f work/qc.csv || gsutil cp "gs://cpg-tob-wgs-${main_suf}/gvcf/batch${batch}/*.csv" work/qc.csv
	R --vanilla <<code
rmarkdown::render('qc.Rmd', output_file='qc.html', params=list(\
gender_tsv='work/gender.tsv', \
age_csv='work/age.csv', \
qc_csv='work/qc.csv', \
gvcf_bucket_suffix='${main_suf}'\
))
code
	gsutil cp qc.html gs://cpg-tob-wgs-${web_suf}/qc/qc.html
}

# Run test first
run 0 test test temporary

# Run on full data in a standard access level
if [[ $PROD = "YES" ]]
then
	run $BATCH analysis main web
fi
