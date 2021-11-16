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
    --new-batches)
	NEWBATCHES="$2"
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
		r-reactable r-ggrepel r-sessioninfo r-gargle r-here r-assertthat r-dt r-googlecloudstorager r-ggforce r-plotly
fi

function run() {
	local namespace=$1
	local joint_calling_run_version=$2
	local new_batches=$3
	local dir=work

	if [[ $TMP = "YES" ]]
	then
		META_TSV_PATH="gs://cpg-tob-wgs-${namespace}-tmp/analysis/joint-calling/${joint_calling_run_version}/meta.tsv"
	else
		META_TSV_PATH="gs://cpg-tob-wgs-${namespace}-analysis/joint-calling/${joint_calling_run_version}/meta.tsv"
	fi
	test -f ${dir}/meta.tsv || gsutil cp ${META_TSV_PATH} ${dir}/meta.tsv

	cat qc.Rmd | sed 's/ fig.width=[1-9]*, fig.height=[1-9]*/ fig.width=plot_width, fig.height=plot_height/g' > qc-for-html.Rmd
	R --vanilla <<code
rmarkdown::render('qc-for-html.Rmd', output_file='qc.html', params=list(\
meta_tsv='${dir}/meta.tsv', \
namespace='${namespace}', \
new_batches='${new_batches}'
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
	run main $VERSION $NEWBATCHES
else
	run test $VERSION $NEWBATCHES
fi

set +x
