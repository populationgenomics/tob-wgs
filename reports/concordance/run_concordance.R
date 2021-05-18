#!/usr/bin/env Rscript

require(tidyverse)
require(glue)

gcs_outdir <- Sys.getenv("OUTPUT")
container_output <- "mtcars.tsv"
x <- mtcars
readr::write_tsv(x, file = container_output)
system("gcloud -q auth activate-service-account --key-file=/gsa-key/key.json")
system(glue("gcloud cp {container_output} {gcs_outdir}"))
head(x)
