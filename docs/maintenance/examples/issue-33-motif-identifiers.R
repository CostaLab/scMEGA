suppressPackageStartupMessages({library(Signac);library(SeuratObject);library(SummarizedExperiment);library(dplyr)})
env <- new.env(parent=globalenv()); env$assay <- SummarizedExperiment::assay
sys.source('R/trajectory_analysis.R',env)
ids <- paste0('M',1:5)
labels <- c('GATA1','GATA1::TAL1','GATA1','Spi1','Actb')
names(labels) <- ids
peaks <- c('chr1-100-200','chr1-300-400','chr2-500-600')
match.matrix <- Matrix::Matrix(matrix(c(1,0,1,0,1,0,1,1,0,0,1,1,1,0,1),3,
 dimnames=list(peaks,ids)),sparse=TRUE)
pwm <- setNames(rep(list(matrix(0.25,4,4,dimnames=list(c('A','C','G','T'),NULL))),length(ids)),ids)
motifs <- Signac::CreateMotifObject(data=match.matrix,pwm=pwm,motif.names=as.list(labels))
counts <- Matrix::Matrix(matrix(1:24,3,dimnames=list(peaks,paste0('cell',1:8))),sparse=TRUE)
chromatin <- Signac::CreateChromatinAssay(counts,motifs=motifs,min.cells=0,min.features=0)
actual.labels <- unlist(Signac::GetMotifData(chromatin,slot='motif.names'))
stopifnot(identical(unname(actual.labels),make.unique(unname(labels))))
cat('PASS: Signac makes duplicate labels unique; suffixed labels are not gene identifiers.\n')
set.seed(33)
activity <- matrix(rnorm(40),5,dimnames=list(ids,paste0('T',1:8)))
expression <- matrix(rnorm(32),4,dimnames=list(c('Gata1','Tal1','Spi1','Actb'),colnames(activity)))
expression['Spi1',] <- activity['M4',]+seq_len(8)/10
expression['Actb',] <- activity['M5',]+seq_len(8)/20
# This reproduces the exact label assignment/matching used by scMEGA, not a species mapper.
rownames(activity) <- unname(actual.labels)
mm <- SummarizedExperiment::SummarizedExperiment(assays=list(mat=activity))
rna <- SummarizedExperiment::SummarizedExperiment(assays=list(mat=expression))
result <- env$GetCorrelation(mm,rna)
stopifnot(setequal(result$tfs,c('Spi1','Actb')),!'Gata1' %in% result$tfs)
cat('PASS: exact matching retains Spi1/Actb but loses mouse Gata1 despite GATA1 motif labels.\n')
case.changed <- stringr::str_to_title(actual.labels)
stopifnot(anyDuplicated(stringr::str_to_title(c('GATA1','Gata1')))>0L,'Gata1::Tal1' %in% case.changed,
          !'Gata1::Tal1' %in% rownames(expression))
cat('PASS: title casing leaves a composite unmatched and can collapse distinct-case GATA1/Gata1 labels.\n')
reordered.ids <- rev(ids)
stopifnot(!identical(unname(actual.labels[reordered.ids]),unname(actual.labels)))
cat('PASS: positional motif-label assignment is unsafe after motif row reordering; stable IDs must be checked.\n')
annotations <- new.env()
load('data/geneAnnoMm10.rda',annotations)
genes <- annotations$geneAnnoMm10$genes
cat('Bundled mm10 symbols:',head(S4Vectors::mcols(genes)$symbol),'\n')
cat('Sequence names:',head(unique(as.character(GenomicRanges::seqnames(genes)))),'\n')
sys.source('R/select_tf_gene.R',env)
stopifnot(!'genome' %in% names(formals(env$SelectTFs)), identical(formals(env$SelectGenes)$genome, 'hg38'))
cat('PASS: SelectTFs has no genome argument; SelectGenes defaults to hg38 and needs explicit mm10.\n')
