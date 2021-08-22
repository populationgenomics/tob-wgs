#!/usr/bin/env bash

set -ex

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c bioconda -c conda-forge plink

# Copy Plink files from main-upload to main,
# and transfer a subset to test.

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

PLINK_UPLOAD_DIR="gs://cpg-tob-wgs-main-upload/one1k1_genotyping_all_regions"
SUBSET_DIR="subset_20samples"
OUTPUT_BUCKET="gs://cpg-tob-wgs-main/snpchip/v1/plink"
gsutil -m cp "${PLINK_UPLOAD_DIR}/*" .
mkdir -p ${SUBSET_DIR}

plink \
    --bfile onek1k \
    --keep ${SAMPLES_FILE} \
    --make-bed \
    --out  ${SUBSET_DIR}/onek1k

gsutil -m cp "${SUBSET_DIR}/*" ${OUTPUT_BUCKET}/

DESCRIPTION="OneK1K SNPchip data in Plink format, provided 2021-Jul-20 by GWCCG,\nsubset to 20 samples for test bucket."
printf ${DESCRIPTION} | gsutil cp - ${OUTPUT_BUCKET}/README.txt
