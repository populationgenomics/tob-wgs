#!/usr/bin/env bash

set -ex

gsutil -m cp cpg-tob-wgs-main/scrna-seq/CellRegMap_input_files/expression_objects/sce22.h5ad gs://cpg-tob-wgs-test/scrna-seq/CellRegMap_input_files/expression_objects/
gsutil -m cp cpg-tob-wgs-main/scrna-seq/CellRegMap_input_files/all_B_cells/PCs_Bcells.csv.pkl gs://cpg-tob-wgs-test/scrna-seq/CellRegMap_input_files/expression_PCs/
