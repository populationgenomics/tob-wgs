#!/usr/bin/env bash

set -ex

gsutil -m cp gs://cpg-tob-wgs-main/CellRegMap_input_files/plink_files/plink_chr22.bed \
    gs://cpg-tob-wgs-main/CellRegMap_input_files/plink_files/plink_chr22.bim \
    gs://cpg-tob-wgs-main/CellRegMap_input_files/plink_files/plink_chr22.fam \
    gs://cpg-tob-wgs-main/CellRegMap_input_files/GRM/grm_wide.csv \
    gs://cpg-tob-wgs-main/CellRegMap_input_files/expression_objects/sce22.h5ad gs://cpg-tob-wgs-test/v0/skat/
