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

# copy VCF + spreadsheet into container
#system(glue("gsutil cp {snpchipid2tobid_excel} {snpchip_vcf_raw} ."))

# get snpchip sample names from VCF

cat("Testing micromamba")
system("micromamba install --name base -c bioconda -c conda-forge bcftools")
system("bcftools --help")
cat(glue("[{date()}] Finished!!!"))
