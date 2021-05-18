#!/usr/bin/env Rscript

require(tidyverse)
require(glue)

# authorise gcp service account
gcp_auth_sa <- "gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"
system(gcp_auth_sa)

gcs_outdir <- Sys.getenv("OUTPUT")
gcs_indir <- "gs://cpg-tob-wgs-snpchipdata/data"
snpchipid2tobid_excel <- glue("{gcs_indir}/OneK1K_sample_IDs_2021-Apr-15.xlsx")
snpchip_vcf_raw <- glue("{gcs_indir}/onek1k_pre_imputation_genotypes.vcf.gz")

system(glue("gsutil cp {snpchipid2tobid_excel} {snpchip_vcf_raw} ."))
cat("Listing files in current working directory:")
list.files()
cat(glue("[{date()}] Finished!!!"))
