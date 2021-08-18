#!/usr/bin/env bash

set -ex

INPUT_PATH="gs://cpg-tob-wgs-main-upload/one1k1_genotyping_all_regions"
OUTPUT_PATH="snpchip/v1/plink"
DESCRIPTION="OneK1K SNPchip data in Plink format, provided 2021-Jul-20 by GWCCG."

# Move to 'main'.
gsutil -m mv ${INPUT_PATH} gs://cpg-tob-wgs-main/${OUTPUT_PATH}
echo ${DESCRIPTION} | gsutil cp - gs://cpg-tob-wgs-main/${OUTPUT_PATH}/README.txt
