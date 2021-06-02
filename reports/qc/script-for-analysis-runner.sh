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
	test -f ${dir}/meta.tsv   || gsutil cp "gs://cpg-tob-wgs-${main_suf}/joint-calling/${joint_calling_run_version}/meta.tsv" ${dir}/meta.csv
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
	qc_fpath=qc/qc-${joint_calling_run_version}.html
	gsutil cp qc.html gs://cpg-tob-wgs-${main_suf}-web/${qc_fpath}

	echo ""
	echo "Copied to the ${main_suf}-web bucket. The report will be available under https://${main_suf}-web.populationgenomics.org.au/tob-wgs/${qc_fpath}"
}

# Run on full data in a standard access level
if [[ $PROD = "YES" ]]
then
	run $BATCH main $VERSION
else
	run batch1 test $VERSION
fi
