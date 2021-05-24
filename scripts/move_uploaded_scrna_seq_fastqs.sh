#!/bin/bash

set -ex

gsutil mv gs://cpg-tob-wgs-upload/H7HTJDSXX gs://cpg-tob-wgs-main/scrna_seq/fastq/v0
gsutil mv gs://cpg-tob-wgs-upload/onek1k_scrna_seq_data_md5sums.txt gs://cpg-tob-wgs-main/scrna_seq/fastq/v0

# Copy a subset of the data to the `test` bucket.
gsutil -m cp -r gs://cpg-tob-wgs-main/scrna_seq/fastq/v0/OneK1K_scRNA_Sample1 gs://cpg-tob-wgs-main/scrna_seq/fastq/v0/OneK1K_scRNA_Sample2 gs://cpg-tob-wgs-test/scrna_seq/fastq/v0

