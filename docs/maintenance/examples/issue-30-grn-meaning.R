suppressPackageStartupMessages(library(dplyr))
e <- new.env(parent=globalenv());sys.source('R/grn.R',e)
m <- Matrix::Matrix(matrix(c(1,1,0,0,1,1),3,dimnames=list(c('p1','p2','p3'),c('TF1','TF2'))),sparse=TRUE)
p2g <- data.frame(peak=c('p1','p2','p3'),gene=c('G1','G1','G2'))
correlations <- data.frame(tf=c('TF1','TF2'),gene=c('G1','G2'),correlation=c(.8,-.6),fdr=c(.01,.02))
r <- e$GetGRN(motif.matching=m,df.cor=correlations,df.p2g=p2g)
stopifnot(r$n_peaks[r$tf=='TF1' & r$gene=='G1']==2,
 r$n_peaks[r$tf=='TF2' & r$gene=='G1']==1,
 is.na(r$correlation[r$tf=='TF2' & r$gene=='G1']),
 r$correlation[r$tf=='TF2' & r$gene=='G2']==-.6)
print(r)
cat('PASS: GetGRN joins supplied evidence, counts motif/peak-link rows, and does not itself remove missing/negative correlations.\n')
cat('No trajectory-free statistical inference or causal validation occurs in this assembly function.\n')
