#!/bin/bash

# Overwrite previous metadata files, as a column had been swapped.
gsutil mv gs://cpg-tob-wgs-upload/R_210315_BINKAN1_1K1KDNA_M001.csv gs://cpg-tob-wgs-main/gvcf/batch0
gsutil mv gs://cpg-tob-wgs-upload/R_210315_BINKAN1_1K1KDNA_M002.csv gs://cpg-tob-wgs-main/gvcf/batch1
