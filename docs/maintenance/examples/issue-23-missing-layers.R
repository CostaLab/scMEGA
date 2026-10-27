suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(SummarizedExperiment);library(dplyr)})
env <- new.env(parent=globalenv())
env$Assays <- SeuratObject::Assays; env$LayerData <- SeuratObject::LayerData
env$SummarizedExperiment <- SummarizedExperiment::SummarizedExperiment
sys.source('R/trajectory_analysis.R',env)
x <- Matrix::Matrix(matrix(1:16,2,dimnames=list(c('g1','g2'),paste0('c',1:8))),sparse=TRUE)
obj <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts=x),assay='SoupXRNA')
obj$Trajectory <- c(10,20,30,40,60,70,80,90)
missing <- tryCatch(env$GetTrajectory(obj,assay='SoupXRNA',slot='data',groupEvery=50,
 scaleTo=NULL,log2Norm=FALSE,smoothWindow=NULL),error=identity)
stopifnot(inherits(missing,'error'))
cat('EXPECTED missing data failure:',conditionMessage(missing),'\n')
obj <- Seurat::NormalizeData(obj,assay='SoupXRNA',verbose=FALSE)
r <- env$GetTrajectory(obj,assay='SoupXRNA',slot='data',groupEvery=50,
 scaleTo=NULL,log2Norm=FALSE,smoothWindow=NULL)
stopifnot(identical(dim(r),c(2L,2L)),all(is.finite(SummarizedExperiment::assay(r))))
cat('PASS: intended RNA normalization provides the required data layer, with no assay casting.\n')
activity <- matrix(seq(-8,7),2,dimnames=list(c('M1','M2'),colnames(obj)))
obj[['chromvar']] <- SeuratObject::CreateAssayObject(data=activity)
r <- env$GetTrajectory(obj,assay='chromvar',slot='data',groupEvery=50,
 scaleTo=NULL,log2Norm=FALSE,smoothWindow=NULL)
expected <- cbind(rowMeans(activity[,1:4,drop=FALSE]),rowMeans(activity[,5:8,drop=FALSE]))
stopifnot(isTRUE(all.equal(unname(SummarizedExperiment::assay(r)),unname(expected))))
cat('PASS: a data-only ordinary Assay can hold TF activity; ChromatinAssay casting is unnecessary.\n')
