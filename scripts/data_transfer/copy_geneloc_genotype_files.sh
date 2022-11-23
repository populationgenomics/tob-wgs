#!/usr/bin/env bash

set -ex

# copy chromosome 7 and 22 geneloc files from main to test
gsutil cp gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr7.tsv gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/
gsutil cp gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/
