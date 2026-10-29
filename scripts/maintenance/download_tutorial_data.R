#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("usage: download_tutorial_data.R <tutorial> <destination-root>")

manifest <- read.delim("scripts/maintenance/tutorial-data.tsv",
                       stringsAsFactors = FALSE, check.names = FALSE)
manifest <- manifest[manifest$tutorial == args[[1]], , drop = FALSE]
if (!nrow(manifest))
  stop("unknown tutorial: ", args[[1]])

options(timeout = 7200)
for (i in seq_len(nrow(manifest))) {
  destination <- file.path(args[[2]], manifest$destination[[i]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(destination) || file.info(destination)$size != manifest$bytes[[i]])
    download.file(manifest$url[[i]], destination, mode = "wb", quiet = FALSE)
  actual <- unname(tools::sha256sum(destination))
  if (!identical(actual, manifest$sha256[[i]]))
    stop("checksum mismatch: ", destination)
  message("verified ", destination)
}
