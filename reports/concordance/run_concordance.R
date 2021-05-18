#!/usr/bin/env Rscript

require(tidyverse)
require(glue)

gcs_outdir <- Sys.getenv("OUTPUT")
mtcars_tsv <- "mtcars.tsv"
mtcars_plot <- "mtcars_disp_vs_hp.png"
d <- mtcars
readr::write_tsv(d, file = container_output)

p <- d %>%
  ggplot(aes(x = disp, y = hp)) +
  geom_point()
ggsave(mtcars_plot, p)

system("gcloud -q auth activate-service-account --key-file=/gsa-key/key.json")
system(glue("gsutil cp {mtcars_tsv} {mtcars_plot} {gcs_outdir}"))
cat("That's all folks!!!")
