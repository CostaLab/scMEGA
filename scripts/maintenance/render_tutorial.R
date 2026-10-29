#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("usage: render_tutorial.R <tutorial> <data-root> <output-dir>")

articles <- c(
  myofibroblast = "vignettes/articles/myofibroblast-GRN.Rmd",
  cardiomyocyte = "vignettes/articles/cardiomyocyte-GRN.Rmd",
  cd4 = "vignettes/articles/pbmc_10x_multiome.Rmd"
)
if (!args[[1]] %in% names(articles))
  stop("unknown tutorial: ", args[[1]])

dir.create(args[[3]], recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  articles[[args[[1]]],
  output_dir = normalizePath(args[[3]], mustWork = TRUE),
  knit_root_dir = normalizePath(args[[2]], mustWork = TRUE),
  envir = new.env(parent = globalenv()),
  clean = TRUE
)
