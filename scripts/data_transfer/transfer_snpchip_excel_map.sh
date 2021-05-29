#!/usr/bin/env bash

# Transfer SNPchip ID-to-TOB ID Excel file into main bucket
INPUT_EXCEL="gs://cpg-tob-wgs-snpchipdata/data/OneK1K_sample_IDs_2021-Apr-15.xlsx"

if [[ "${OUTPUT}" =~ "gs://cpg-tob-wgs-main/snpchip/v1" ]]; then
    gsutil cp "${INPUT_EXCEL}" "${OUTPUT}"
fi
