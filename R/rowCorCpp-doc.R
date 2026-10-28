#' Correlate selected pairs of matrix rows
#'
#' Compute Pearson correlations between corresponding rows selected from two
#' dense numeric matrices. This is the compiled helper used by PeakToGene.
#'
#' @name rowCorCpp
#' @usage rowCorCpp(idxX, idxY, X, Y)
#' @param idxX Non-empty integer vector of one-based row indices into X.
#' @param idxY Integer vector of one-based row indices into Y, with the same length
#' as idxX. Each element selects the row paired with the corresponding idxX row.
#' @param X Dense numeric matrix with observations in columns.
#' @param Y Dense numeric matrix with the same observations in the same column
#' order as X. The function checks column counts, not observation identifiers.
#' @details Callers must provide valid positive in-range indices of equal length.
#' Correlations for constant rows can be non-finite. The helper does not perform
#' normalization, imputation, multiple-testing adjustment, or identifier matching.
#' @return A numeric vector with one correlation for each pair of row indices.
#' @seealso \code{\link{PeakToGene}}
NULL
