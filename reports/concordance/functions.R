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

save_obj_rds <- function(o, b, dir, date) {
  saveRDS(o, glue("{dir}/{b}/{date}_list_contents.rds"))
}

read_obj_rds <- function(b, dir, date) {
  readRDS(glue("{dir}/{b}/{date}_list_contents.rds"))
}

chipid2tobid <- function(x, id_map) {
  stopifnot(length(x) == 1, c("sample_id", "tob_id") %in% names(id_map))
  tobid <- id_map %>%
    dplyr::filter(sample_id == x)
  if (nrow(tobid) == 1) {
    tobid$tob_id
  } else if (nrow(tobid) == 0) {
    message("No matches! Returning NA.")
    return(NA_character_)
  } else {
    stop("Wait, how did you get more than 1 matches??")
  }
}

tobid2chipid <- function(x, id_map) {
  stopifnot(length(x) == 1, c("sample_id", "tob_id") %in% names(id_map))
  chipid <- id_map %>%
    dplyr::filter(tob_id == x)
  if (nrow(chipid) == 1) {
    chipid$sample_id
  } else if (nrow(chipid) == 0) {
    message("No matches! Returning NA.")
    return(NA_character_)
  } else {
    stop("Wait, how did you get more than 1 matches??")
  }
}

run_gvcftools <- function(input, output) {
  idir <- normalizePath(dirname(input))
  ibname <- basename(input)
  dir.create(dirname(output), recursive = TRUE, showWarnings = TRUE)
  odir <- normalizePath(dirname(output))
  obname <- basename(output)

  docker_cmd <- glue(
    "docker run --rm ",
    "-v {idir}:/data/input -v {odir}:/data/output ",
    "quay.io/biocontainers/gvcftools:0.17.0--he941832_3 /bin/bash -c ")
  gvcftools_cmd <- glue("'gzip -dc /data/input/{ibname} | extract_variants | gzip -c > /data/output/{obname}'")

  system(glue("{docker_cmd} {gvcftools_cmd}"))
}

concordance_stats <- function(wgs, snp) {
  gt_match <- snp %>%
    left_join(wgs, by = "chrpos") %>%
    mutate(gt_equal = gt.x == gt.y)

    # filter(chrpos %in% snp[["chrpos"]]) %>%
    # mutate(gt = sub("\\|", "/", gt)) %>%
    # left_join(snp, by = "chrpos") %>%
    # mutate(gt_match = gt.x == gt.y)

  nrow_snp <- nrow(snp)
  nrow_wgs <- nrow(wgs)
  n_wgs_in_snp <- sum(wgs[["chrpos"]] %in% snp[["chrpos"]])
  n_snp_in_wgs <- sum(snp[["chrpos"]] %in% wgs[["chrpos"]])
  pct_snp_in_wgs <- n_snp_in_wgs / nrow_snp
  pct_wgs_in_snp <- n_wgs_in_snp / nrow_wgs
  n_gt_match <- gt_match %>% filter(gt_equal == TRUE) %>% nrow()
  n_gt_nomatch <- gt_match %>% filter(gt_equal == FALSE) %>% nrow()
  pct_gt_match <- n_gt_match / nrow_snp
  pct_gt_nomatch <- n_gt_nomatch / nrow_snp

  list(
    summary = tibble::tibble(
      pct_gt_match = pct_gt_match,
      pct_gt_nomatch = pct_gt_nomatch,
      nrow_snp = nrow_snp,
      nrow_wgs = nrow_wgs,
      n_wgs_in_snp = n_wgs_in_snp,
      n_snp_in_wgs = n_snp_in_wgs,
      pct_snp_in_wgs = pct_snp_in_wgs,
      pct_wgs_in_snp = pct_wgs_in_snp,
      n_gt_match = n_gt_match,
      n_gt_nomatch = n_gt_nomatch
    ),
    gt_match = gt_match
  )
}

is_outlier <- function(x, m = 1.5) {
  (x < quantile(x, 0.25) - m * IQR(x)) |
    (x > quantile(x, 0.75) + m * IQR(x))
}

select_snpchip_sample <- function(snpchip, tobid, id_map) {
  snpchip %>%
    select(chrpos, ref, alt, gt = tobid2chipid(tobid, id_map)) %>%
    filter(!gt %in% "0/0") # homref
}

read_wgs_vcf <- function(wgsvcf, tobid) {
  data.table::fread(
    cmd = glue("zgrep -v '^##' {wgsvcf}"), sep = "\t", na.strings = ".",
    drop = c("QUAL", "ID", "FILTER", "INFO", "FORMAT"), data.table = FALSE
  ) %>%
    rename(chr = `#CHROM`, pos = POS, ref = REF, alt = ALT) %>%
    mutate(
      gt = sub("(^.*?):.*", "\\1", .data[[tobid]]),
      chrpos = glue("{chr}_{pos}")) %>%
    select(chrpos, ref, alt, gt) %>%
    tibble::as_tibble()
}

count_gt_per_chrom <- function(d, chrom) {
  # Instead of processing full dataset,
  # split by chromosome (non-parallel).
  d %>%
    filter(chr %in% chrom) %>%
    tidyr::pivot_longer(-chr) %>%
    dplyr::group_by(chr, name) %>%
    dplyr::count(value)
}
