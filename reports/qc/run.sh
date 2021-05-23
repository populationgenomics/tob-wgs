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

function run() {
	local batch=$1
	local analysis_suf=$2
	local main_suf=$3
	local web_suf=$4

	test -f qc.csv || gsutil cp "gs://cpg-tob-wgs-main/gvcf/batch${batch}/*.csv" qc.csv
	test -f gender.tsv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/gender.tsv gender.tsv
	test -f age.csv || gsutil cp gs://cpg-tob-wgs-${analysis_suf}/age.csv age.csv
	Rscript -e "rmarkdown::render('qc.Rmd', output_file='qc.html', params=list(gender_tsv='gender.tsv', age_csv='age.csv', qc_csv='qc.csv', bucket_suffix='${main_suf}'))"
	gsutil cp qc.html gs://cpg-tob-wgs-${web_suf}/qc/qc.html
}

# Run test first
run 0 test test temporary

# Run on full data in a standard access level
if [ $PROD -eq YES ]
then
	run $BATCH analysis main web
fi
