#!/usr/bin/env bash

set -ex

gsutil cp "gs://cpg-tob-wgs-main/concordance/v1/v[1-2]-raw_chr22_samples.tsv" "gs://cpg-tob-wgs-test/$OUTPUT"
