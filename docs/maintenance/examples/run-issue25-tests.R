# Focused source tests; package namespace checks are recorded separately.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1]])
suppressPackageStartupMessages({library(testthat); library(SeuratObject); library(dplyr); library(SummarizedExperiment)})
env <- new.env(parent = globalenv())
# Match scMEGA's intended imports; SummarizedExperiment also exports an Assays name.
env$Assays <- SeuratObject::Assays
env$LayerData <- SeuratObject::LayerData
env$SimpleList <- S4Vectors::SimpleList
env$SummarizedExperiment <- SummarizedExperiment::SummarizedExperiment
env$assay <- SummarizedExperiment::assay
sys.source(file.path(root, "R", "utils.R"), env)
sys.source(file.path(root, "R", "trajectory_analysis.R"), env)
results <- testthat::test_dir(file.path(root,"tests","testthat"), filter="trajectory-dimensions",
  env=env, load_package="none", reporter="summary", stop_on_failure=FALSE)
s <- as.data.frame(results)
print(s[,c("test","failed","error","skipped","warning","passed")],row.names=FALSE)
if(any(s$failed > 0 | s$error | s$skipped)) stop("Trajectory regression failure/skip")
