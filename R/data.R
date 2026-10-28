#' Bundled hg19 gene annotation
#'
#' Historical annotation used by \code{PeakToGene} when \code{genome = "hg19"}.
#' These objects are distributed with scMEGA; this documentation does not update
#' their identifiers, coordinates, or annotation release.
#'
#' @format A \code{GenomicRanges::GRangesList} with three elements:
#' \describe{
#'   \item{genes}{Gene ranges with \code{gene_id} and \code{symbol} metadata.}
#'   \item{exons}{Exon ranges.}
#'   \item{TSS}{Transcription start site ranges.}
#' }
#' @source The annotation distributed in the scMEGA repository at baseline commit
#' \code{e6ce4f9}. The objects do not contain an annotation release/accession;
#' the historical preparation provenance requires separate verification.
#' @seealso \code{\link{PeakToGene}}
#' @keywords datasets
"geneAnnoHg19"

#' Bundled hg38 gene annotation
#'
#' Historical annotation used by \code{PeakToGene} when \code{genome = "hg38"}.
#' These objects are distributed with scMEGA; this documentation does not update
#' their identifiers, coordinates, or annotation release.
#'
#' @format A \code{GenomicRanges::GRangesList} with three elements:
#' \describe{
#'   \item{genes}{Gene ranges with \code{gene_id} and \code{symbol} metadata.}
#'   \item{exons}{Exon ranges.}
#'   \item{TSS}{Transcription start site ranges.}
#' }
#' @source The annotation distributed in the scMEGA repository at baseline commit
#' \code{e6ce4f9}. The objects do not contain an annotation release/accession;
#' the historical preparation provenance requires separate verification.
#' @seealso \code{\link{PeakToGene}}
#' @keywords datasets
"geneAnnoHg38"

#' Bundled mm9 gene annotation
#'
#' Historical annotation used by \code{PeakToGene} when \code{genome = "mm9"}.
#' These objects are distributed with scMEGA; this documentation does not update
#' their identifiers, coordinates, or annotation release.
#'
#' @format A \code{GenomicRanges::GRangesList} with three elements:
#' \describe{
#'   \item{genes}{Gene ranges with \code{gene_id} and \code{symbol} metadata.}
#'   \item{exons}{Exon ranges.}
#'   \item{TSS}{Transcription start site ranges.}
#' }
#' @source The annotation distributed in the scMEGA repository at baseline commit
#' \code{e6ce4f9}. The objects do not contain an annotation release/accession;
#' the historical preparation provenance requires separate verification.
#' @seealso \code{\link{PeakToGene}}
#' @keywords datasets
"geneAnnoMm9"

#' Bundled mm10 gene annotation
#'
#' Historical annotation used by \code{PeakToGene} when \code{genome = "mm10"}.
#' These objects are distributed with scMEGA; this documentation does not update
#' their identifiers, coordinates, or annotation release.
#'
#' @format A \code{GenomicRanges::GRangesList} with three elements:
#' \describe{
#'   \item{genes}{Gene ranges with \code{gene_id} and \code{symbol} metadata.}
#'   \item{exons}{Exon ranges.}
#'   \item{TSS}{Transcription start site ranges.}
#' }
#' @source The annotation distributed in the scMEGA repository at baseline commit
#' \code{e6ce4f9}. The objects do not contain an annotation release/accession;
#' the historical preparation provenance requires separate verification.
#' @seealso \code{\link{PeakToGene}}
#' @keywords datasets
"geneAnnoMm10"
