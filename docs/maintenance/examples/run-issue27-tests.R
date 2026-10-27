args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[[1]])
suppressPackageStartupMessages({library(testthat);library(Seurat);library(Signac);library(dplyr)})
env <- new.env(parent=globalenv())
env$LayerData <- SeuratObject::LayerData
env$JoinLayers <- SeuratObject::JoinLayers
sys.source(file.path(root,'R','pair_cells.R'),env)
r <- testthat::test_dir(file.path(root,'tests','testthat'),filter='coembed-layers',env=env,
  load_package='none',reporter='summary',stop_on_failure=FALSE)
s <- as.data.frame(r)
print(s[,c('test','failed','error','skipped','warning','passed')],row.names=FALSE)
if(any(s$failed>0|s$error|s$skipped))stop('Coembedding regression failure/skip')
