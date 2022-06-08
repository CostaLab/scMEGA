SelectTFs <- function(object,
                      motif.assay = "chromvar",
                      rna.assay = "RNA",
                      trajectory.name = "Trajectory",
                      p.cutoff = 0.01,
                      cor.cutoff = 0.3
                      ){
  trajMM <- GetTrajectory(object,
                          assay = motif.assay,
                          slot = "data",
                          smoothWindow = 7,
                          log2Norm = FALSE)

  rownames(trajMM) <- object@assays$ATAC@motifs@motif.names

  trajRNA <- GetTrajectory(object,
                           assay = rna.assay,
                           slot = "data",
                           smoothWindow = 7,
                           log2Norm = TRUE)

  df.cor <- GetCorrelation(trajMM, trajRNA)

  # we only select TFs that show significant correlation
  df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > cor.cutoff, ]

  matMM <- TrajectoryHeatmap(trajMM,
                             varCutOff = 0,
                             pal = paletteContinuous(set = "solarExtra"),
                             limits = c(-2, 2),
                             name = "TF activity",
                             returnMatrix = TRUE)

  df_tf_time_point <- data.frame(tfs = rownames(matMM),
                                 time_point = seq(1, 100, length.out = nrow(matMM)))
  rownames(df_tf_time_point) <- df_tf_time_point$tfs

  df_tf_time_point <- df_tf_time_point[df.cor$tfs, ]
  df.cor$time_point <- df_tf_time_point$time_point
  df.cor <- df.cor[order(df.cor$time_point), ]

  trajMM <- trajMM[df.cor$tfs, ]
  trajRNA <- trajRNA[df.cor$tfs, ]

  ht <- CorrelationHeatmap(trajectory1 = trajMM,
                           trajectory2 = trajRNA,
                           name1 = "TF activity",
                           name2 = "Gene expression")


}
