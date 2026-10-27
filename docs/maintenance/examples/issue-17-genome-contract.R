suppressPackageStartupMessages({library(GenomicRanges);library(SummarizedExperiment);library(S4Vectors)})
e <- new.env(parent=globalenv())
for (name in c('Hg19','Hg38','Mm9','Mm10')) load(paste0('data/geneAnno',name,'.rda'),e)
sys.source('R/peak_to_gene.R',e)
sys.source('R/select_tf_gene.R',e)
Rcpp::sourceCpp('src/Correlation.cpp',env=e)
ann <- e$geneAnnoMm10$genes
symbols <- S4Vectors::mcols(ann)$symbol
unique.symbols <- names(which(table(symbols)==1L))
genes <- ann[symbols %in% unique.symbols][1:3]
symbols <- S4Vectors::mcols(genes)$symbol
tss <- ifelse(as.character(strand(genes))=='+',start(genes),end(genes))
peak.names <- paste0(as.character(seqnames(genes)),'-',tss+101L,'-',tss+201L)
set.seed(17)
gene.mat <- matrix(runif(24,1,20),3,dimnames=list(symbols,paste0('T',1:8)))
peak.mat <- gene.mat + matrix(runif(24,0,.1),3)
rownames(peak.mat) <- peak.names
gene.mat <- gene.mat[c(3,1,2),,drop=FALSE]
unsupported <- tryCatch(e$PeakToGene(peak.mat,gene.mat,genome='rn6'),error=identity)
stopifnot(inherits(unsupported,'error'),grepl('hg19, hg38, mm9, and mm10',conditionMessage(unsupported)))
cat('EXPECTED unsupported genome:',conditionMessage(unsupported),'\n')
r <- e$PeakToGene(peak.mat,gene.mat,genome='mm10',max.dist=5000)
stopifnot(nrow(r)>=3L,setequal(r$gene,symbols),all(r$peak %in% peak.names),all(r$distance>=0))
for(i in seq_len(nrow(r))) {
 expected <- cor(peak.mat[r$peak[i],],gene.mat[r$gene[i],])
 stopifnot(abs(r$Correlation[i]-expected)<1e-10)
}
cat('PASS:',nrow(r),'mm10 candidate links retain peak/gene identities and match direct correlations after gene-row permutation.\n')
stopifnot(identical(formals(e$SelectGenes)$genome,'hg38'), identical(formals(e$PeakToGene)$genome,'hg19'))
cat('PASS: SelectGenes and direct PeakToGene have different genome defaults; explicit assembly selection is required.\n')
