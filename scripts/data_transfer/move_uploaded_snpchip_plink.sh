#!/usr/bin/env bash

# Copy Plink files from main-upload to main,
# and transfer a subset of 20 samples to test.

set -ex

# install Plink
micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c bioconda -c conda-forge plink

SAMPLES_FILE="samples_to_keep.txt"

cat > ${SAMPLES_FILE} << EOF
39 39
40 40
41 41
42 42
43 43
44 44
45 45
46 46
47 47
48 48
49 49
50 50
51 51
52 52
53 53
54 54
55 55
56 56
57 57
58 58
EOF

INPUT_PATH="gs://cpg-tob-wgs-main-upload/one1k1_genotyping_all_regions"
OUTPUT_PATH="snpchip/v1/plink"
OUTPUT_TEST="gs://cpg-tob-wgs-test/${OUTPUT_PATH}"
OUTPUT_MAIN="gs://cpg-tob-wgs-main/${OUTPUT_PATH}"
SUBSET_DIR="subset_20samples"

gsutil -m cp "${INPUT_PATH}/*" .
mkdir -p ${SUBSET_DIR}

plink \
    --bfile onek1k \
    --keep ${SAMPLES_FILE} \
    --make-bed \
    --out  ${SUBSET_DIR}/onek1k

gsutil -m cp "${SUBSET_DIR}/*" ${OUTPUT_TEST}

README_TEST="OneK1K SNPchip data in Plink format, provided 2021-Jul-20 by GWCCG (subset to 20 samples for test bucket)."
echo ${README_TEST} | gsutil cp - ${OUTPUT_TEST}/README.txt

# Move to 'main'.
gsutil -m mv ${INPUT_PATH} ${OUTPUT_MAIN}
README_MAIN="OneK1K SNPchip data in Plink format, provided 2021-Jul-20 by GWCCG."
echo ${README_MAIN} | gsutil cp - ${OUTPUT_MAIN}/README.txt
