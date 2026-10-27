# Real Seurat/Signac pipeline, small synthetic counts, no biological download.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else getwd()
suppressPackageStartupMessages({library(Seurat);library(Signac);library(dplyr)})
env <- new.env(parent=globalenv())
env$LayerData <- SeuratObject::LayerData
env$JoinLayers <- SeuratObject::JoinLayers
sys.source(file.path(root,'R','pair_cells.R'),env)
set.seed(27)
ng <- 100; nr <- 120; nq <- 100
loadings <- matrix(runif(ng*4,0.2,3),ng,4)
rna.latent <- matrix(runif(4*nr,0.2,4),4,nr)
atac.latent <- matrix(runif(4*nq,0.2,4),4,nq)
make.counts <- function(mu,rows,cols) {
 m <- matrix(as.numeric(rpois(length(mu),mu)),nrow(mu),dimnames=list(rows,cols))
 Matrix::Matrix(m,sparse=TRUE)
}
rna <- make.counts(loadings %*% rna.latent,paste0('gene',1:ng),paste0('r',1:nr))
ga <- make.counts(loadings %*% atac.latent,rownames(rna),paste0('a',1:nq))
peaks <- make.counts(matrix(runif(ng*4,.2,3),ng,4)%*%atac.latent,paste0('chr1-',seq(100,10000,100),'-',seq(150,10050,100)),colnames(ga))
query <- CreateSeuratObject(CreateChromatinAssay(peaks,min.cells=0,min.features=0),assay='ATAC')
query <- RunTFIDF(query)
query <- FindTopFeatures(query,min.cutoff=0)
query <- RunSVD(query,n=30)
ref <- CreateSeuratObject(CreateAssay5Object(counts=list(batch1=rna[,1:60],batch2=rna[,61:120])),assay='RNA')
ref <- NormalizeData(ref,verbose=FALSE)
ref <- FindVariableFeatures(ref,nfeatures=ng,verbose=FALSE)
ref.joined <- JoinLayers(ref,assay='RNA')
stopifnot(isTRUE(all.equal(as.matrix(LayerData(ref.joined,assay='RNA',layer='counts'))[,colnames(rna)],as.matrix(rna))))
ref.legacy <- ref.joined
ref.legacy[['RNA']] <- as(ref.legacy[['RNA']], 'Assay')
for (kind in c('split','joined','legacy')) {
 cat('\nCASE:',kind,'\n')
 object <- switch(kind, split=ref, joined=ref.joined, legacy=ref.legacy)
 result <- tryCatch(env$CoembedData(object,query,ga,weight.reduction='lsi',dims=1:10,verbose=FALSE),error=function(e)e)
 if(inherits(result,'error')) {
   cat('ERROR:',conditionMessage(result),'\n')
 } else {
   stopifnot(setequal(colnames(result),c(colnames(rna),colnames(ga))),all(is.finite(Embeddings(result,'umap'))))
   cat('PASS:',ncol(result),'cells with finite UMAP. RNA class:',class(result[['RNA']])[1],'\n')
 }
}
