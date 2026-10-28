# Standalone reproducer: real SeuratObject accessor, no scMEGA installation needed.
suppressPackageStartupMessages(library(SeuratObject))
x <- Matrix::Matrix(matrix(seq_len(24), nrow = 2,
  dimnames = list(c("geneA", "geneB"), paste0("c", seq_len(12)))), sparse = TRUE)
requested <- rev(colnames(x))
a <- CreateAssay5Object(counts = x)
observed <- LayerData(a, layer = "counts", cells = requested)
cat("Requested:", requested, "\nReturned:", colnames(observed), "\n")
original <- observed
colnames(original) <- paste0("pair", seq_len(ncol(original)))
expected <- x[, requested, drop = FALSE]
colnames(expected) <- colnames(original)
stopifnot(!identical(as.matrix(original), as.matrix(expected)))
corrected <- observed[, requested, drop = FALSE]
colnames(corrected) <- colnames(original)
stopifnot(identical(as.matrix(corrected), as.matrix(expected)))
cat("CONFIRMED: original renaming loses cell identity; explicit name indexing preserves counts.\n")
print(sessionInfo())
