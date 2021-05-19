#!/usr/bin/env Rscript

require(dplyr, include.only = c("mutate", "select", "%>%"))
require(glue, include.only = "glue")
require(readr, include.only = "write_tsv")
require(readxl, include.only = "read_excel")
require(stringr, include.only = "str_order")
require(tibble, include.only = "as_tibble_col")

# authorise gcp service account
gcp_auth_sa <- "gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"
system(gcp_auth_sa)

cat("Installing bcftools\n")
system("micromamba install --name base -c bioconda -c conda-forge bcftools")

gcs_outdir <- Sys.getenv("OUTPUT")
gcs_indir <- "gs://cpg-tob-wgs-snpchipdata/data"
snpchipid2tobid_excel <- glue("{gcs_indir}/OneK1K_sample_IDs_2021-Apr-15.xlsx")
snpchip_vcf_raw <- glue("{gcs_indir}/onek1k_pre_imputation_genotypes.vcf.gz")

# copy VCF + spreadsheet into container
system(glue("gsutil cp {snpchipid2tobid_excel} {snpchip_vcf_raw} ."))

# 1,034 samples in SNPchip VCF
chip_sample_nms <-
  system(glue("{bcftools} query -l {snpchip_vcf_raw}"), intern = TRUE) %>%
  as_tibble_col("sample")

# 983 samples in Excel spreadsheet
chipid2tobid <-
  snpchipid2tobid_excel %>%
  read_excel(skip = 1) %>%
  mutate(tob_id_final = glue("TOB{TOB_ID}")) %>%
  select(sample = PERSON, tob_id = tob_id_final)

# TRUE - just write the map with 983 samples to use
# for VCF reheading. The 51 SNPchip samples without
# a TOB ID will stay with their SNPchip IDs.
all(chipid2tobid$sample %in% chip_sample_nms$sample)

chipid2tobid_tsv <- "chipid2tobid.tsv"
# order by TOB ID
ord <- str_order(chipid2tobid$tob_id, numeric = TRUE)
chipid2tobid[ord, ] %>%
  write_tsv(file = chipid2tobid_tsv, col_names = FALSE)

# rename VCF samples to TOB IDs according to map
vcf_rehead <- "0-rehead.vcf.bgz"
system(glue("bcftools reheader ",
            "-s {chipid2tobid_tsv} ",
            "-o {vcf_rehead} ",
            "{snpchip_vcf_raw}"))

# move data to bucket
system(glue("gsutil cp {vcf_rehead} {gcs_outdir}"))

cat(glue("[{date()}] Finished!!!"))
