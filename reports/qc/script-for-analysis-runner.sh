#!/bin/bash

set -ex

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -a|--analysis-runner)
    ANALYSIS_RUNNER=YES
    shift # past argument
    ;;
    -p|--prod)
    PROD=YES
    shift # past argument
    ;;
    --tmp)
    TMP=YES
    shift
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

if [[ $ANALYSIS_RUNNER = "YES" ]] ; then
	micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
		r-reactable r-ggrepel r-sessioninfo r-gargle r-here r-assertthat r-dt r-googlecloudstorager r-ggforce
fi

function run() {
	local batch=$1
	local namespace=$2
	local joint_calling_run_version=$3
	local dir=work

	test -f ${dir}/reported_sex.tsv || gsutil cp "gs://cpg-tob-wgs-${namespace}-metadata/reported_sex.tsv" ${dir}/reported_sex.tsv
	test -f ${dir}/age.csv          || gsutil cp "gs://cpg-tob-wgs-${namespace}-metadata/age.csv" ${dir}/age.csv
	test -f ${dir}/qc.csv           || gsutil cp "gs://cpg-tob-wgs-${namespace}-metadata/${batch}/*.csv" ${dir}/qc.csv
	if [[ $TMP = "YES" ]]
	then
		test -f ${dir}/meta.tsv     || gsutil cp "gs://cpg-tob-wgs-test-tmp/joint-calling/${joint_calling_run_version}/meta.tsv" ${dir}/meta.tsv
	else
		test -f ${dir}/meta.tsv     || gsutil cp "gs://cpg-tob-wgs-${namespace}-metadata/joint-calling/${joint_calling_run_version}/meta.tsv" ${dir}/meta.tsv
	fi
	cat qc.Rmd | sed 's/r fig.width=[1-9]*, fig.height=[1-9]*/r fig.width=plot_width, fig.height=plot_height/g' > qc-for-html.Rmd
	R --vanilla <<code
rmarkdown::render('qc-for-html.Rmd', output_file='qc.html', params=list(\
reported_sex_tsv='${dir}/reported_sex.tsv', \
age_csv='${dir}/age.csv', \
qc_csv='${dir}/qc.csv', \
meta_tsv='${dir}/meta.tsv', \
namespace='${namespace}'
))
code
	rm qc-for-html.Rmd
	qc_fpath=qc/qc-${joint_calling_run_version}.html
	gsutil cp qc.html gs://cpg-tob-wgs-${namespace}-web/${qc_fpath}

	echo ""
	echo "Copied to the ${namespace}-web bucket. The report will be available under https://${namespace}-web.populationgenomics.org.au/tob-wgs/${qc_fpath}"
}

if [[ $PROD = "YES" ]]
then
	run $BATCH main $VERSION
else
	run $BATCH test $VERSION
fi
