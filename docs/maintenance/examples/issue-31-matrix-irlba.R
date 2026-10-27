# Run in a fresh R process; use R_LIBS_USER to choose the tested irlba build.
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(irlba))
cat('R:', as.character(getRversion()), 'Matrix:', as.character(packageVersion('Matrix')),
    'irlba:', as.character(packageVersion('irlba')), '\n')
cat('Matrix path:', find.package('Matrix'), '\nirlba path:', find.package('irlba'), '\n')
set.seed(31)
a <- Matrix::rsparsematrix(80, 40, density=0.2, rand.x=function(n) runif(n,1,10))
x <- irlba::irlba(a, nv=3, tol=1e-8, maxit=1000)
reference <- svd(as.matrix(a), nu=0, nv=0)$d[1:3]
stopifnot(max(abs(x$d-reference)) < 1e-6,
          max(abs(crossprod(x$v)-diag(3))) < 1e-6)
relative.residual <- norm(as.matrix(a %*% x$v)-sweep(x$u,2,x$d,'*'),'F') / norm(as.matrix(a),'F')
stopifnot(relative.residual < 1e-6)
cat('PASS: sparse truncated SVD matches dense reference; residual',relative.residual,'\n')

suppressPackageStartupMessages(library(Signac))
set.seed(31)
counts <- Matrix::rsparsematrix(120,60,density=0.25,rand.x=function(n) as.numeric(sample.int(8,n,replace=TRUE)))
rownames(counts) <- paste0('chr1-', seq(100,12000,100), '-', seq(150,12050,100))
colnames(counts) <- paste0('cell',seq_len(60))
obj <- SeuratObject::CreateSeuratObject(Signac::CreateChromatinAssay(counts, min.cells=0, min.features=0), assay='ATAC')
obj <- Signac::RunTFIDF(obj)
obj <- Signac::FindTopFeatures(obj, min.cutoff=0)
obj <- Signac::RunSVD(obj, n=5)
embedding <- SeuratObject::Embeddings(obj,'lsi')
stopifnot(identical(dim(embedding), c(60L,5L)),all(is.finite(embedding)),identical(rownames(embedding),colnames(counts)))
cat('PASS: Signac TF-IDF / feature selection / SVD preserves cell identities.\n')
