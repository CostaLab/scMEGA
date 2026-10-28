#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("usage: benchmark_coembed_memory.R <input-fixture-rds> <output-summary-rds>")

suppressPackageStartupMessages({
  library(scMEGA)
  library(Seurat)
  library(Signac)
})

set.seed(1234)
fixture <- readRDS(args[[1]])
Rprofmem(paste0(args[[2]], ".Rprofmem"), threshold = 1e6)
started <- proc.time()
result <- CoembedData(
  obj.rna = fixture$rna,
  obj.atac = fixture$atac,
  gene.activity = fixture$activity,
  weight.reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)
elapsed <- proc.time() - started
Rprofmem(NULL)

summary <- list(
  elapsed = unname(elapsed),
  cells = colnames(result),
  metadata = result[[]],
  variable.features = VariableFeatures(result),
  assays = Assays(result),
  geneactivity.scale = LayerData(
    result,
    assay = "GeneActivity",
    layer = "scale.data"
  ),
  rna.data = LayerData(result, assay = "RNA", layer = "data"),
  pca = Embeddings(result, "pca"),
  umap = Embeddings(result, "umap")
)
saveRDS(summary, args[[2]], version = 3)
print(elapsed)
