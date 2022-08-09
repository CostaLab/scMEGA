#' Select TFs for GRN inference
#'
#' @param object A Seurat object
#' @param tf.assay The assay name for TF activity. Default: "chromvar"
#' @param rna.assay The assay name for gene expression. Default: "RNA"
#' @param trajectory.name The trajectory name used for computing correlation
#' between TF binding activity and TF expression
#' @param p.cutoff A cutoff of p-values. Default: 0.01
#' @param cor.cutoff A cutoff of correlation. Default: 0.3
#' @param return.heatmap Whether or not return the heatmap for visualization
#'
#' @return A list containing a dataframe of TF activity and expression correlation
#'and a heatmap
#' @export
#'
SelectTFs <- function(object,
                      tf.assay = "chromvar",
                      rna.assay = "RNA",
                      trajectory.name = "Trajectory",
                      p.cutoff = 0.01,
                      cor.cutoff = 0.3,
                      return.heatmap = TRUE) {
  
    trajMM <- suppressMessages(GetTrajectory(
    object,
    assay = tf.assay,
    slot = "data",
    smoothWindow = 7,
    log2Norm = FALSE
  )) 

  rownames(trajMM) <- object@assays$ATAC@motifs@motif.names

  trajRNA <- suppressMessages(GetTrajectory(
    object,
    assay = rna.assay,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  ))

  df.cor <- GetCorrelation(trajMM, trajRNA)

  # we only select TFs that show significant correlation
  df.cor <-
    df.cor[df.cor$adj_p < p.cutoff &
             df.cor$correlation > cor.cutoff,]

  matMM <- suppressMessages(TrajectoryHeatmap(
    trajMM,
    varCutOff = 0,
    pal = paletteContinuous(set = "solarExtra"),
    limits = c(-2, 2),
    name = "TF activity",
    returnMatrix = TRUE
  ))

  df_tf_time_point <- data.frame(tfs = rownames(matMM),
                                 time_point = seq(1, 100, length.out = nrow(matMM)))
  rownames(df_tf_time_point) <- df_tf_time_point$tfs

  df_tf_time_point <- df_tf_time_point[df.cor$tfs,]
  df.cor$time_point <- df_tf_time_point$time_point
  df.cor <- df.cor[order(df.cor$time_point),]

  trajMM <- trajMM[df.cor$tfs,]
  trajRNA <- trajRNA[df.cor$tfs,]


  if (return.heatmap) {
    ht <- suppressMessages(
      CorrelationHeatmap(
        trajectory1 = trajMM,
        trajectory2 = trajRNA,
        name1 = "TF activity",
        name2 = "Gene expression"
      )
    )

    res <- list("tfs" = df.cor, "heatmap" = ht)

  } else{
    res <- list("tfs" = df.cor)

  }

  return(res)

}

#' Select genes
#'
#' @param object A Seurat object
#' @param atac.assay The assay name for chromatin accessibility. Default: "ATAC"
#' @param rna.assay The assay name for gene expression. Default: "RNA"
#' @param var.cutoff.gene The cutoff of variation to select genes. Default: 0.9
#' @param trajectory.name The trajectory name used for computing correlation
#' between TF binding activity and TF expression
#' @param distance.cutoff The minimum distance between cis-regulatory elements and genes
#' @param cor.cutoff The cutoff of peak-to-gene correlation. Default: 0
#' @param fdr.cutoff The cutoff of peak-to-gene p-value Default: 1e-04
#' @param return.heatmap Whether or not return the heatmap for visualization
#' @param labelTop1 Number of labels for row names
#' @param labelTop2 Number of labels for row names
#'
#' @return A list containing a dataframe of peak-to-gene links and a heatmap
#' @export
#'
SelectGenes <- function(object,
                        atac.assay = "ATAC",
                        rna.assay = "RNA",
                        var.cutoff.gene = 0.9,
                        trajectory.name = "Trajectory",
                        distance.cutoff = 2000,
                        cor.cutoff = 0,
                        fdr.cutoff = 1e-04,
                        return.heatmap = TRUE,
                        labelTop1 = 10,
                        labelTop2 = 10) {
  trajRNA <- GetTrajectory(
    object,
    assay = rna.assay,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  )

  trajATAC <- GetTrajectory(
    object,
    assay = atac.assay,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  )

  # note here we only use the top 10% most variable genes
  groupMatRNA <- suppressMessages(
    TrajectoryHeatmap(
      trajRNA,
      varCutOff = var.cutoff.gene,
      pal = paletteContinuous(set = "horizonExtra"),
      limits = c(-2, 2),
      returnMatrix = TRUE
    )
  )

  groupMatATAC <- suppressMessages(
    TrajectoryHeatmap(
      trajATAC,
      varCutOff = 0,
      maxFeatures = nrow(trajATAC),
      pal = paletteContinuous(set = "solarExtra"),
      limits = c(-2, 2),
      name = "Chromatin accessibility",
      returnMatrix = TRUE
    )
  )

  message("Linking cis-regulatory elements to genes...")
  df.p2g <- PeakToGene(peak.mat = groupMatATAC,
                       gene.mat = groupMatRNA,
                       genome = "hg38")

  df.p2g <- df.p2g %>%
    subset(distance > distance.cutoff) %>%
    subset(Correlation > cor.cutoff & FDR < fdr.cutoff)

  trajATAC <- trajATAC[df.p2g$peak,]
  trajRNA <- trajRNA[df.p2g$gene,]

  if (return.heatmap) {
    ht <- suppressMessages(
      CorrelationHeatmap(
        trajectory1 = trajATAC,
        trajectory2 = trajRNA,
        name1 = "Chromatin accessibility",
        name2 = "Gene expression",
        labelTop1 = labelTop1,
        labelTop2 = labelTop2,
        labelRows1 = FALSE,
        labelRows2 = FALSE
      )
    )

    res <- list("p2g" = df.p2g, "heatmap" = ht)

  } else{
    res <- list("p2g" = df.p2g)
  }


  return(res)

}
