guess_file_type <- function(x) {
  dplyr::case_when(
    grepl("\\.bam$", x, ignore.case = TRUE) ~ "BAM",
    grepl("\\.bai$", x, ignore.case = TRUE) ~ "BAMindex",
    grepl("\\.cram$", x, ignore.case = TRUE) ~ "CRAM",
    grepl("\\.crai$", x, ignore.case = TRUE) ~ "CRAMindex",
    grepl("\\.fastq.gz$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("\\.fastq$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("\\.fq$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("\\.fq\\.gz$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("manifest\\.txt$", x, ignore.case = TRUE) ~ "Manifest",
    grepl("\\.md5$", x, ignore.case = TRUE) ~ "MD5",
    grepl("md5\\.txt$", x, ignore.case = TRUE) ~ "MD5txt",
    grepl("\\.vcf$", x, ignore.case = TRUE) ~ "VCF_unz",
    grepl("\\.g\\.vcf\\.gz$", x, ignore.case = TRUE) ~ "GVCF",
    grepl("\\.vcf\\.gz$", x, ignore.case = TRUE) ~ "VCF",
    grepl("\\.tbi$", x, ignore.case = TRUE) ~ "VCFindex",
    grepl("\\.csv$", x, ignore.case = TRUE) ~ "CSV",
    TRUE ~ "OTHER")
}

my_gcs_list_obj <- function(b) {
  googleCloudStorageR::gcs_list_objects(bucket = b,
                                        detail = "summary") %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(name = glue("gs://{b}/{name}"),
                  size = sub(" bytes", "", size), # else returns NA
                  size = fs::as_fs_bytes(size),
                  ftype = guess_file_type(name))
}

save_obj_rds <- function(o, b) {
  saveRDS(o, glue("{dir_bucket}/{b}/{date}_list_contents.rds"))
}

read_obj_rds <- function(b) {
  readRDS(glue("{dir_bucket}/{b}/{date}_list_contents.rds"))
}
