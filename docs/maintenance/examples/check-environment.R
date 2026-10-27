# Detect drift before re-running the recorded source validations; no installation.
args <- commandArgs(trailingOnly=TRUE)
manifest <- if(length(args))args[1] else 'docs/maintenance/phase2-evidence/environment-manifest.tsv'
x <- read.delim(manifest,stringsAsFactors=FALSE,check.names=FALSE,na.strings='')
stopifnot(as.character(getRversion())=='4.5.3')
actual <- installed.packages()
versions <- actual[match(x$Package,rownames(actual)),'Version']
bad <- is.na(versions) | versions != x$Version
if(any(bad)) {
 print(data.frame(package=x$Package[bad],expected=x$Version[bad],actual=versions[bad]))
 stop('Validation environment changed; reproduce/record the new stack before claiming the same validation.')
}
for(k in seq_len(nrow(x))) {
 root <- file.path(actual[match(x$Package[k],rownames(actual)),'LibPath'],x$Package[k])
 stopifnot(unname(tools::md5sum(file.path(root,'DESCRIPTION')))==x$DESCRIPTION_MD5[k])
 files <- sort(list.files(file.path(root,'libs'),pattern='[.](so|dll|dylib)$',full.names=TRUE,recursive=TRUE))
 checksum <- paste(paste(basename(files),unname(tools::md5sum(files)),sep=':'),collapse=';')
 expected <- if(is.na(x$BINARY_MD5[k])) '' else x$BINARY_MD5[k]
 stopifnot(identical(checksum,expected))
}
cat('PASS: R and',nrow(x),'recorded package versions, descriptions and shared-library checksums match.\n')
