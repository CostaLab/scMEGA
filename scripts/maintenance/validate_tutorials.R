#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork = TRUE)
articles <- c(
  myofibroblast = "vignettes/articles/myofibroblast-GRN.Rmd",
  cardiomyocyte = "vignettes/articles/cardiomyocyte-GRN.Rmd",
  cd4 = "vignettes/articles/pbmc_10x_multiome.Rmd"
)

required <- list(
  myofibroblast = c("CoembedData", "CreatePairedObject", "SelectTFs",
                    "SelectGenes", "GetTFGeneCorrelation", "GetGRN"),
  cardiomyocyte = c("CoembedData", "CreatePairedObject", "SelectTFs",
                    "SelectGenes", "GetTFGeneCorrelation", "GetGRN"),
  cd4 = c("SelectTFs", "SelectGenes", "GetTFGeneCorrelation", "GetGRN")
)

namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))

for (name in names(articles)) {
  path <- file.path(root, articles[[name]])
  stopifnot(file.exists(path))
  source <- readLines(path, warn = FALSE)
  if (any(grepl("--no-check-certificate", source, fixed = TRUE)))
    stop(name, " disables TLS verification")
  code <- tempfile(fileext = ".R")
  knitr::purl(path, output = code, quiet = TRUE)
  parse(code)
  missing <- setdiff(required[[name]], exports)
  if (length(missing))
    stop(name, " references unexported scMEGA functions: ", paste(missing, collapse = ", "))
  for (fun in required[[name]]) {
    if (!any(grepl(paste0("\\b", fun, "\\s*\\("), source)))
      stop(name, " no longer exercises ", fun, "()")
  }
  message(name, ": parse and API checks passed")
}

cardiomyocyte <- paste(
  readLines(file.path(root, articles[["cardiomyocyte"]]), warn = FALSE),
  collapse = "\n"
)
sampling.contract <- c(
  "USE_DOWNSAMPLING <- TRUE",
  "N_RNA_CELLS <- 40000L",
  "N_ATAC_CELLS <- 8000L",
  "SAMPLING_SEED <- 42L",
  "patient_group",
  "patient_region_id",
  "gene.activity <- gene.activity[, selected.atac.cells, drop = FALSE]"
)
missing.contract <- sampling.contract[!vapply(
  sampling.contract,
  grepl,
  logical(1),
  x = cardiomyocyte,
  fixed = TRUE
)]
if (length(missing.contract))
  stop("Cardiomyocyte sampling contract is incomplete: ",
       paste(missing.contract, collapse = ", "))
message("Cardiomyocyte fixed sampling contract passed")

pkgdown <- paste(readLines(file.path(root, "_pkgdown.yml"), warn = FALSE), collapse = "\n")
for (article in basename(articles)) {
  href <- sub("\\.Rmd$", ".html", article)
  if (!grepl(href, pkgdown, fixed = TRUE))
    stop("_pkgdown.yml does not link ", href)
}

manifest <- read.delim(file.path(root, "scripts/maintenance/tutorial-data.tsv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(manifest) == 10L)
stopifnot(all(manifest$tutorial %in% names(articles)))
stopifnot(all(grepl("^[0-9a-f]{64}$", manifest$sha256)))
stopifnot(all(grepl("^https://", manifest$url)))
stopifnot(!anyDuplicated(manifest$destination))
message("Tutorial manifest and pkgdown navigation checks passed")
