#!/usr/bin/env Rscript

if (!requireNamespace("BPCells", quietly = TRUE)) {
  stop("Install optional BPCells before running this evaluation.")
}
suppressPackageStartupMessages({
  library(scMEGA)
  library(Seurat)
  library(Signac)
})

set.seed(71)
root <- tempfile("scmega-bpcells-")
dir.create(root)
counts <- Matrix::rsparsematrix(100, 40, density = 0.08)
counts@x <- as.numeric(abs(round(counts@x * 5)))
counts <- methods::as(counts, "dgCMatrix")
rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
BPCells::write_matrix_dir(
  BPCells::convert_matrix_type(counts, "uint32_t"),
  file.path(root, "rna")
)
bp <- BPCells::open_matrix_dir(file.path(root, "rna"))
stopifnot(identical(dimnames(bp), dimnames(counts)))
stopifnot(all(as.matrix(bp) == as.matrix(counts)))

bp.object <- SeuratObject::CreateSeuratObject(
  SeuratObject::CreateAssay5Object(counts = bp),
  assay = "RNA"
)
mem.object <- SeuratObject::CreateSeuratObject(
  SeuratObject::CreateAssay5Object(counts = counts),
  assay = "RNA"
)
bp.object <- Seurat::NormalizeData(bp.object, verbose = FALSE)
mem.object <- Seurat::NormalizeData(mem.object, verbose = FALSE)
normalized.diff <- max(abs(
  as.matrix(SeuratObject::LayerData(bp.object, layer = "data")) -
    as.matrix(SeuratObject::LayerData(mem.object, layer = "data"))
))

trajectory <- seq(2.5, 100, length.out = ncol(bp.object))
bp.object$Trajectory <- trajectory
mem.object$Trajectory <- trajectory
bp.traj <- GetTrajectory(
  bp.object,
  assay = "RNA",
  slot = "data",
  groupEvery = 20,
  smoothWindow = NULL,
  log2Norm = FALSE
)
mem.traj <- GetTrajectory(
  mem.object,
  assay = "RNA",
  slot = "data",
  groupEvery = 20,
  smoothWindow = NULL,
  log2Norm = FALSE
)
trajectory.diff <- max(abs(
  SummarizedExperiment::assay(bp.traj) -
    SummarizedExperiment::assay(mem.traj)
))

legacy.result <- try(
  SeuratObject::CreateAssayObject(counts = bp),
  silent = TRUE
)
peaks <- counts[seq_len(50), ]
rownames(peaks) <- paste0(
  "chr1-", seq(100, 5000, 100), "-", seq(150, 5050, 100)
)
BPCells::write_matrix_dir(
  BPCells::convert_matrix_type(peaks, "uint32_t"),
  file.path(root, "atac")
)
bp.peaks <- BPCells::open_matrix_dir(file.path(root, "atac"))
chromatin <- Signac::CreateChromatinAssay(
  counts = bp.peaks,
  min.cells = 0,
  min.features = 0
)
cat("BPCells", as.character(utils::packageVersion("BPCells")), "\n")
cat("raw_counts_max_abs 0\n")
cat("normalized_max_abs", normalized.diff, "\n")
cat("trajectory_max_abs", trajectory.diff, "\n")
cat("legacy_CreateAssayObject_supported", !inherits(legacy.result, "try-error"), "\n")
if (!inherits(legacy.result, "try-error")) {
  cat("legacy_counts_class",
      paste(class(SeuratObject::LayerData(legacy.result, layer = "counts")), collapse = "/"),
      "\n")
}
cat("counts_layer_class",
    paste(class(SeuratObject::LayerData(bp.object, layer = "counts")), collapse = "/"),
    "\n")
cat("chromatin_counts_class",
    paste(class(SeuratObject::LayerData(chromatin, layer = "counts")), collapse = "/"),
    "\n")
