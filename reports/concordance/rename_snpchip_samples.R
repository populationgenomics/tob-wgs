#!/usr/bin/env Rscript

# Rename SNPchip VCF sample names based on a TSV map.
#
# Input:
#  - VCF with SNPchip genotypes
#  - Excel file with mapping between SNPchip ID and TOB ID
# Output:
#  - VCF with samples renamed

require(dplyr, include.only = c("mutate", "select", "%>%"))
require(glue, include.only = "glue")
require(readr, include.only = "write_tsv")
require(readxl, include.only = "read_excel")
require(stringr, include.only = "str_order")
require(tibble, include.only = "as_tibble_col")

cat(glue("[{date()}] Installing bcftools.\n"))
system("micromamba install --name base -c bioconda -c conda-forge bcftools")

# set up inputs
gcs_indir <- "gs://cpg-tob-wgs-main/snpchip/v1"
snpchip_vcf_raw_gcs <- glue("{gcs_indir}/onek1k_pre_imputation_genotypes.vcf.gz")
snpchip_vcf_raw <- basename(snpchip_vcf_raw_gcs)
snpchipid2tobid_excel_gcs <- glue("{gcs_indir}/OneK1K_sample_IDs_2021-Apr-15.xlsx")
snpchipid2tobid_excel <- basename(snpchipid2tobid_excel_gcs)

# copy inputs to container
system(glue("gsutil cp {snpchipid2tobid_excel_gcs} {snpchip_vcf_raw_gcs} ."))

# 1,034 samples in SNPchip VCF
chip_sample_nms <-
  system(glue("bcftools query -l {snpchip_vcf_raw}"), intern = TRUE) %>%
  as_tibble_col("sample")

# 983 samples in Excel file
chipid2tobid <-
  snpchipid2tobid_excel %>%
  read_excel(skip = 1) %>%
  mutate(tob_id = glue("TOB{TOB_ID}")) %>%
  select(sample = PERSON, tob_id)

# all SNPchip IDs in the Excel file correspond to SNPchip IDs
# in the VCF - just write the map with 983 samples to use
# for VCF reheading. The 51 SNPchip samples without
# a TOB ID will stay with their SNPchip IDs.
all(chipid2tobid$sample %in% chip_sample_nms$sample)

chipid2tobid_tsv <- "chipid2tobid.tsv"
# order by TOB ID
ord <- str_order(chipid2tobid$tob_id, numeric = TRUE)
chipid2tobid[ord, ] %>%
  write_tsv(file = chipid2tobid_tsv, col_names = FALSE)

# rename VCF samples to TOB IDs according to map
snpchip_vcf_rehead <- "snpchip_vcf_rehead.vcf.bgz"
system(glue("bcftools reheader ",
            "-s {chipid2tobid_tsv} ",
            "-o {snpchip_vcf_rehead} ",
            "{snpchip_vcf_raw}"))

# copy output to analysis bucket
gcs_outdir <- Sys.getenv("OUTPUT")
stopifnot(grepl("gs://cpg-tob-wgs-analysis/snpchip", gcs_outdir))
system(glue("gsutil cp {snpchip_vcf_rehead} {gcs_outdir}"))

cat(glue("[{date()}] Finished!!!\n"))
