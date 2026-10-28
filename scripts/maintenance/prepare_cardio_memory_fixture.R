#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("usage: prepare_cardio_memory_fixture.R <cardiomyocyte-data-dir> <output-rds>")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

rna <- readRDS(file.path(args[[1]], "snRNA.rds"))
atac <- readRDS(file.path(args[[1]], "snATAC.rds"))
activity <- readRDS(file.path(args[[1]], "gene.activity.rds"))
genes <- sort(intersect(rownames(activity), rownames(rna)))
genes <- genes[seq_len(min(2000L, length(genes)))]
rna.cells <- sort(colnames(rna))[seq_len(1200L)]
atac.cells <- sort(colnames(atac))[seq_len(600L)]

fixture <- list(
  rna = subset(rna, features = genes, cells = rna.cells),
  atac = subset(atac, cells = atac.cells),
  activity = activity[genes, atac.cells, drop = FALSE],
  genes = genes,
  rna.cells = rna.cells,
  atac.cells = atac.cells
)
saveRDS(fixture, args[[2]], version = 3)
