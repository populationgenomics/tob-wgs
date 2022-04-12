#!/usr/bin/env bash

set -ex

gsutil cp "gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/covariates_files/B_*_peer_factors_file.txt" gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/