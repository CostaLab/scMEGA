#' Get peak to gene links
#'
#' This function will link peak to genes based on correlation analysis.
#'
#' @param peak.mat A matrix containing chromatin accessibility data
#' @param gene.mat A matrix containing gene expression data
#' @param genome Which genome to use. Currently available are: hg19, hg38, mm9,
#' and mm10
#' @param max.dist The maximum distance between a peak and a gene
#' @useDynLib scMEGA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @return A data frame containing all peak-to-gene links
#' @export
#'
PeakToGene <- function(peak.mat,
                       gene.mat,
                       genome = "hg19",
                       max.dist = 250000) {
  if (!(genome %in% c("hg19", "hg38", "mm9", "mm10"))) {
    stop("Available genome are: hg19, hg38, mm9, and mm10!")
  }

  if (genome == "hg19") {
    gene_anno <- geneAnnoHg19
  } else if (genome == "hg38") {
    gene_anno <- geneAnnoHg38
  } else if (genome == "mm9") {
    gene_anno <- geneAnnoMm9
  } else if (genome == "mm10") {
    gene_anno <- geneAnnoMm10
  }

  ## create object for RNA data
  genes <- gene_anno$genes
  gene.use <-
    intersect(elementMetadata(genes)[, "symbol"], rownames(gene.mat))

  genes <- genes[elementMetadata(genes)[, "symbol"] %in% gene.use]
  gene.mat <- gene.mat[gene.use,]

  gene_start <- ifelse(genes@strand == "+",
                       genes@ranges@start,
                       genes@ranges@start + genes@ranges@width - 1)

  genes <- GRanges(
    genes@seqnames,
    ranges = IRanges(gene_start,
                     width = 1),
    name = genes$symbol,
    gene_id = genes$gene_id,
    strand = genes@strand
  )

  seRNA <-
    SummarizedExperiment(assays = SimpleList(RNA = gene.mat),
                         rowRanges = genes)

  ## create object for ATAC data
  df_peak <-
    stringr::str_split_fixed(rownames(groupMatATAC), "-", 3)

  peakSet <- GRanges(df_peak[, 1],
                     IRanges(start = as.numeric(df_peak[, 2]),
                             end = as.numeric(df_peak[, 3])))

  seATAC <-
    SummarizedExperiment(assays = SimpleList(ATAC = peak.mat),
                         rowRanges = peakSet)

  ## find putative peak-to-gene
  o <-
    data.frame(findOverlaps(
      resize(seRNA, 2 * max.dist + 1, "center"),
      resize(rowRanges(seATAC), 1, "center"),
      ignore.strand = TRUE
    ))
  o$distance <- IRanges::distance(rowRanges(seRNA)[o[, 1]],
                                  rowRanges(seATAC)[o[, 2]])
  colnames(o) <- c("gene_idx", "peak_idx", "distance")

  df <- rowRanges(seATAC)[o$peak_idx,]

  o$gene <- rowData(seRNA)[o$gene_idx,]$name
  o$peak <- paste0(
    df@seqnames,
    "-",
    as.data.frame(df@ranges)$start,
    "-",
    as.data.frame(df@ranges)$end
  )


  ## compute correlation
  o$Correlation <- rowCorCpp(as.integer(o$peak_idx),
                             as.integer(o$gene_idx),
                             assay(seATAC),
                             assay(seRNA))

  ## compute p-value
  o$TStat <-
    (o$Correlation / sqrt((
      pmax(1 - o$Correlation ^ 2, 0.00000000000000001, na.rm = TRUE)
    ) / (ncol(seATAC) - 2))) #T-statistic P-value

  o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o <- o[!is.na(o$FDR),]

  return(o)

}
