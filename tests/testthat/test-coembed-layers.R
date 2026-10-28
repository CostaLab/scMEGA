coembed_fixture <- function() {
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
  query <- SeuratObject::CreateSeuratObject(Signac::CreateChromatinAssay(peaks,min.cells=0,min.features=0),assay='ATAC')
  query <- Signac::RunTFIDF(query)
  query <- Signac::FindTopFeatures(query,min.cutoff=0)
  query <- Signac::RunSVD(query,n=30)
  ref <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts=list(batch1=rna[,1:60],batch2=rna[,61:120])),assay='RNA')
  ref <- Seurat::NormalizeData(ref,verbose=FALSE)
  ref <- Seurat::FindVariableFeatures(ref,nfeatures=60,verbose=FALSE)
  ref.joined <- SeuratObject::JoinLayers(ref,assay='RNA')
  ref.scaled <- Seurat::ScaleData(ref.joined,verbose=FALSE)
  stopifnot(isTRUE(all.equal(as.matrix(SeuratObject::LayerData(ref.joined,assay='RNA',layer='counts'))[,colnames(rna)],as.matrix(rna))))
  ref.legacy <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssayObject(counts=rna),assay='RNA')
  ref.legacy <- Seurat::NormalizeData(ref.legacy,verbose=FALSE)
  ref.legacy <- Seurat::FindVariableFeatures(ref.legacy,nfeatures=60,verbose=FALSE)
  list(rna=rna, peaks=peaks, ga=ga, query=query, split=ref,
    joined=ref.joined, scaled=ref.scaled, legacy=ref.legacy)
}

test_that("split reference layers fail with preparation instructions", {
  f <- coembed_fixture()
  expect_error(CoembedData(f$split,f$query,f$ga,weight.reduction="lsi",dims=1:10,verbose=FALSE),
    "one unambiguous 'data' layer")
})

test_that("joined Assay5 and legacy Assay coembeddings retain cells and observed data", {
  f <- coembed_fixture()
  for (kind in c("joined","scaled","legacy")) {
    reference <- f[[kind]]
    expected <- SeuratObject::LayerData(reference,assay="RNA",layer="data")
    coembed.warnings <- character()
    result <- withCallingHandlers(
      CoembedData(reference,f$query,f$ga,weight.reduction="lsi",dims=1:10,verbose=FALSE),
      warning = function(w) {
        coembed.warnings <<- c(coembed.warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expected.patterns <- c(
      "Layer counts isn't present in the assay object",
      "computing too large a percentage of total singular values",
      "default method for RunUMAP has changed",
      "Different features in new layer data than already exists for scale.data",
      "Different cells in new layer data than already exists for scale.data"
    )
    expect_true(any(grepl(expected.patterns[[1]], coembed.warnings, fixed = TRUE)))
    expect_true(
      all(vapply(
        coembed.warnings,
        function(message) any(vapply(
          expected.patterns,
          grepl,
          logical(1),
          x = message,
          fixed = TRUE
        )),
        logical(1)
      )),
      info = paste(coembed.warnings, collapse = " | ")
    )
    expect_setequal(colnames(result),c(colnames(f$rna),colnames(f$ga)))
    expect_equal(ncol(result),220L)
    expect_identical(SeuratObject::VariableFeatures(result), SeuratObject::VariableFeatures(reference[["RNA"]]))
    expect_true(all(is.finite(SeuratObject::Embeddings(result,"umap"))))
    expect_identical(rownames(SeuratObject::Embeddings(result,"umap")),colnames(result))
    actual <- SeuratObject::LayerData(result,assay="RNA",layer="data")
    expect_equal(as.matrix(actual)[rownames(expected),colnames(expected),drop=FALSE],as.matrix(expected))
    expect_equal(as.matrix(SeuratObject::LayerData(result,assay="ATAC",layer="counts"))[,colnames(f$peaks),drop=FALSE],as.matrix(f$peaks))
    expect_identical(rownames(actual),rownames(f$rna))
  }
})
