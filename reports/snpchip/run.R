#!/usr/bin/env Rscript

# CLI for snpchip_ids.Rmd

r_pkgs <- paste0("r-", c("argparser", "dt", "janitor", "glue", "sessioninfo"))
#cat(paste0("[",  date(), "] ", "Installing conda pkgs for report.\n"))
#system(paste("micromamba install --name base -c bioconda -c conda-forge", "bcftools", r_pkgs))

spsm <- suppressPackageStartupMessages
spsm(require(argparser, include.only = c("parse_args", "arg_parser", "add_argument")))
spsm(require(rmarkdown, include.only = "render"))
spsm(require(glue, include.only = "glue"))

p <- arg_parser(
  description = "Exploration of SNPchip and TOB IDs.",
  name = "snpchip", hide.opts = TRUE)
p <- add_argument(p,
                  arg = "--namespace",
                  help = "test or main namespace.")
p <- add_argument(p,
                  arg = "--rmd",
                  help = "Path to input Rmd.")
p <- add_argument(p,
                  arg = "--html",
                  help = "Path to output HTML.")
a <- parse_args(p)

params_list <- list(
  namespace = a$namespace
)

cat("params_list:\n")
print(params_list)

cat(glue("[{as.character(Sys.time())}] START analysis!"))

render(
  input = a$rmd,
  output_file = I(a$html),
  params = params_list
)

cat(glue("[{as.character(Sys.time())}] END concordance!"))
