#!/usr/bin/env bash

set -ex

gsutil -m mv gs://cpg-tob-wgs-main-upload/open_chromatin_annotation/predicted_l1_celltypes_avg_peaks.csv gs://cpg-tob-wgs-main/tob_wgs_rv/open_chromatin_annotation/

gsutil -m mv gs://cpg-tob-wgs-main-upload/open_chromatin_annotation/predicted_l1_celltypes_avg_peaks_chr21.csv gs://cpg-tob-wgs-test/tob_wgs_rv/open_chromatin_annotation/
