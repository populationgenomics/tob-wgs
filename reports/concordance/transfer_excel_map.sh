#!/usr/bin/env bash

# Transfer SNPchip ID-to-TOB ID Excel file into temporary bucket
# gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"

INPUT_EXCEL="gs://cpg-tob-wgs-snpchipdata/data/OneK1K_sample_IDs_2021-Apr-15.xlsx"

if [[ "${OUTPUT}" =~ "gs://cpg-tob-wgs-temporary/peterd-snpchip" ]]; then
    gsutil cp "${INPUT_EXCEL}" "${OUTPUT}"
fi
