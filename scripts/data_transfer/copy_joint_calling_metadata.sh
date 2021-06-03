#!/usr/bin/env bash

VERSION=$1

if [[ -z "$VERSION" ]]; then
    echo "Provide a version tag of a joint-calling run as a first argument" 1>&2
    exit 1
fi

target_bucket=gs://cpg-tob-wgs-test-metadata/joint-calling
if [[ ${OUTPUT} != ${target_bucket}* ]]; then
	echo "The analysis runner output directory must start with ${target_bucket}"
	exit 1
fi	

gsutil cp \
	gs://cpg-tob-wgs-test/joint-calling/${VERSION}/meta.tsv \
	${target_bucket}/${VERSION}/meta.tsv
