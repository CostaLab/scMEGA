SelectTFs <- function(object,
                      motif.assay = "chromvar",
                      rna.assay = "RNA",
                      trajectory.name = "Trajectory",
                      p.cutoff = 0.01,
                      cor.cutoff = 0.3,
                      return.heatmap = TRUE
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

    
if(return.heatmap){
      ht <- suppressMessages(CorrelationHeatmap(trajectory1 = trajMM,
                           trajectory2 = trajRNA,
                           name1 = "TF activity",
                           name2 = "Gene expression"))
    
    res <- list("tfs" = df.cor, "heatmap" = ht)

} else{
    res <- list("tfs" = df.cor)

}

    return(res)

}

SelectGenes <- function(object,
                        atac.assay = "ATAC",
                      rna.assay = "RNA",
                        var.cutoff.gene = 0.9,
                      trajectory.name = "Trajectory",
                       distance.cutoff = 2000,
                       cor.cutoff = 0,
                       fdr.cutoff = 1e-04,
                       return.heatmap = TRUE){
    
    trajRNA <- GetTrajectory(object, 
                        assay = rna.assay, 
                        slot = "data",
                       smoothWindow = 7,
                       log2Norm = TRUE)
    
    trajATAC <- GetTrajectory(object, 
                        assay = atac.assay, 
                        slot = "data",
                       smoothWindow = 7,
                       log2Norm = TRUE)
    
    # note here we only use the top 10% most variable genes
    groupMatRNA <- TrajectoryHeatmap(trajRNA,
                            varCutOff = var.cutoff.gene,
                            pal = paletteContinuous(set = "horizonExtra"),
                            limits = c(-2, 2),
                                     returnMatrix = TRUE)
    
    groupMatATAC <- TrajectoryHeatmap(trajATAC, 
                                      varCutOff = 0,
                                  maxFeatures = nrow(trajATAC),
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                       name = "Chromatin accessibility",
                                      returnMatrix = TRUE)
    
    print(head(groupMatATAC))
    
    df.p2g <- PeakToGene(peak.mat = groupMatATAC,
                     gene.mat = groupMatRNA, 
                     genome = "hg38")
    
    df.p2g <- df.p2g %>%
    subset(distance > distance.cutoff) %>%
    subset(Correlation > cor.cutoff & FDR < fdr.cutoff)
    
    trajATAC <- trajATAC[df.p2g$peak, ]
    trajRNA <- trajRNA[df.p2g$gene, ]

    if(return.heatmap){
            ht <- suppressMessages(CorrelationHeatmap(trajectory1 = trajATAC, 
                             trajectory2 = trajRNA,
                            name1 = "Chromatin accessibility",
                            name2 = "Gene expression",
                            labelRows1 = FALSE,
                            labelRows2 = FALSE))
        
         res <- list("p2g" = df.p2g, "heatmap" = ht)
    
    } else{
        res <- list("p2g" = df.p2g)
    }
    

    return(res)

}
