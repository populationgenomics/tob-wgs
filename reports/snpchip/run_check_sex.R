#!/usr/bin/env Rscript

sort(unname(installed.packages()[,1]))

r_pkgs <- paste("r-dt", "r-sessioninfo", collapse=" ")
cat(paste0("[",  date(), "] ", "Installing conda pkgs for report.\n"))
system(paste("micromamba install --name base -c bioconda -c conda-forge", "r-base=4.0.3", r_pkgs))

spsm <- suppressPackageStartupMessages
spsm(require(rmarkdown, include.only = "render"))
spsm(require(glue, include.only = "glue"))
spsm(require(dplyr))
spsm(require(assertthat, include.only = "assert_that"))

bucket_path <- function(path, bucket_category = NULL) {

  dataset <- Sys.getenv("DATASET")
  access_level <- Sys.getenv("ACCESS_LEVEL")
  assert_that(nchar(dataset) > 0, nchar(access_level) > 0)

  namespace <- dplyr::if_else(access_level == "test", "test", "main")
  if (is.null(bucket_category)) {
    bucket_category <- namespace
  } else if (!bucket_category %in% c("archive", "upload")) {
    bucket_category <- glue("{namespace}-{bucket_category}")
  }

  return(file.path("gs:/", glue("cpg-{dataset}-{bucket_category}"), path))
}

output_path <- function(path_suffix, bucket_category = NULL) {
  output <- Sys.getenv("OUTPUT")
  assert_that(nchar(output) > 0)
  return(bucket_path(file.path(output, path_suffix), bucket_category))
}

params_list <- list(
  input_plink_fam = bucket_path("snpchip/v1/plink/onek1k.fam"),
  input_sex_meta = bucket_path("reported_sex.tsv", "metadata")
)
html_suffix <- "check_sex1.html"
html <- output_path(html_suffix, "web")

cat("params_list:\n")
print(params_list)
cat("html:\n")
print(html)

cat(glue("[{as.character(Sys.time())}] START analysis!"))
render(
  input = "check_sex.Rmd",
  output_file = html_suffix,
  params = params_list
)
system(glue("gsutil cp {html_suffix} {html}"))
cat(glue("[{as.character(Sys.time())}] END analysis!"))
