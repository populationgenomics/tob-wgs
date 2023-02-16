#!/usr/bin/env bash

set -ex

gcloud storage --billing-project tob-wgs cp --recursive \
    gs://cpg-tob-wgs-main/MRR_cram_extracts/2023_01_30 \
    gs://cpg-tob-wgs-release/MRR_cram_extracts/2023_01_30
