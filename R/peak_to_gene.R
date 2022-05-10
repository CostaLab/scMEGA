PeakToGene <- function(peak.mat, gene.mat, genome = "hg19",
                      max.dist = 250000){
    if(!(genome %in% c("hg19", "hg38", "mm9", "mm10"))){
        stop("Available genome are: hg19, hg38, mm9, and mm10!")
    }
    
    if(genome == "hg19"){
        gene_anno <- geneAnnoHg19
    } else if(genome == "hg38"){
        gene_anno <- geneAnnoHg38
    } else if(genome == "mm9"){
        gene_anno <- geneAnnoMm9
    } else if(genome == "mm10"){
        gene_anno <- geneAnnoMm10
    }

    groupMatRNA <- plotTrajectoryHeatmap(trajGEX,
                        varCutOff = 0.9,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2),
                           returnMatrix = TRUE)
    
    seRNA <- SummarizedExperiment(assays = SimpleList(RNA = groupMatRNA), 
                              rowRanges = geneStart)
    

}