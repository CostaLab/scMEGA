# Record the exact dependency closure used by the focused validations.
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args))args[1] else 'docs/maintenance/phase2-evidence'
seed <- c('ArchR','Seurat','SeuratObject','Signac','Matrix','irlba','testthat','SummarizedExperiment',
 'GenomicRanges','S4Vectors','IRanges','dplyr','tidyr','ggraph','igraph','graphlayouts',
 'ComplexHeatmap','destiny','ggpubr','factoextra','optmatch','pracma','nabor','textshape')
i <- installed.packages(fields=c('Repository','RemoteType','RemoteRepo','RemoteSha'))
deps <- tools::package_dependencies(seed,db=i,which=c('Depends','Imports','LinkingTo'),recursive=TRUE)
p <- sort(unique(c(seed,unlist(deps))))
p <- intersect(p,rownames(i))
fields <- c('Package','Version','LibPath','Built','Repository','RemoteType','RemoteRepo','RemoteSha')
x <- as.data.frame(i[match(p,rownames(i)),fields,drop=FALSE],stringsAsFactors=FALSE)
x$DESCRIPTION_MD5 <- vapply(seq_len(nrow(x)),function(k) unname(tools::md5sum(file.path(x$LibPath[k],x$Package[k],'DESCRIPTION'))),character(1))
x$BINARY_MD5 <- vapply(seq_len(nrow(x)),function(k) {
 files <- sort(list.files(file.path(x$LibPath[k],x$Package[k],'libs'),pattern='[.](so|dll|dylib)$',full.names=TRUE,recursive=TRUE))
 paste(paste(basename(files),unname(tools::md5sum(files)),sep=':'),collapse=';')
},character(1))
write.table(x,file.path(out,'environment-manifest.tsv'),quote=TRUE,sep='\t',row.names=FALSE,na='')
writeLines(c(paste('R:',R.version.string),paste('Platform:',R.version$platform),
 paste('Library paths:',paste(.libPaths(),collapse='; ')),
 'This records a tested source-validation environment; it is not a certified full scMEGA/Bioconductor support matrix.'),file.path(out,'environment.txt'))
cat('Recorded',nrow(x),'packages in the validation dependency closure.\n')
