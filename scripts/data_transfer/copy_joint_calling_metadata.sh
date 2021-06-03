#!/usr/bin/env bash

VERSION=$1

if [[ -z "$VERSION" ]]; then
    echo "Provide a version tag of a joint-calling run as a first argument" 1>&2
    exit 1
fi

gsutil cp \
	gs://cpg-tob-wgs-main/joint-calling/${VERSION}/meta.tsv \
	gs://cpg-tob-wgs-main-metadata/joint-calling/${VERSION}/meta.tsv
