GetTFGeneCor <- function(object, 
                         tf.use = NULL,
                         gene.use = NULL,
                         tf.assay = "chromvar", 
                         gene.assay = "RNA",
                        trajectory.name = "Trajectory"){
    
    ## get tf activity and gene expression along trajectory
    trajMM <- GetTrajectory(object, 
                        assay = tf.assay, 
                        slot = "data",
                            trajectory.name = trajectory.name,
                       smoothWindow = 7,
                       log2Norm = FALSE)

    trajRNA <- GetTrajectory(object, 
                        assay = gene.assay, 
                        slot = "data",
                             trajectory.name = trajectory.name,
                       smoothWindow = 7,
                       log2Norm = TRUE)
    
    rownames(trajMM) <- object@assays$ATAC@motifs@motif.names

    tf_activity <- suppressMessages(TrajectoryHeatmap(trajMM, 
                            varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                           name = "TF activity",
                           returnMatrix = TRUE)) 

    gene_expression <- suppressMessages(TrajectoryHeatmap(trajRNA, 
                            varCutOff = 0.9,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                           name = "Gene expression",
                                        returnMatrix = TRUE))
    
    ## here we filter the TFs according to our correlation analysis
    if(!is.null(tf.use)){
        tf_activity <- tf_activity[tf.use, ]
    }
    
    ## here we further filter the genes by only consider genes that are linked to peaks
    if(!is.null(gene.use)){
        sel_genes <- intersect(rownames(gene_expression), gene.use)
        
        ## subset the gene expression matrix
        gene_expression <- gene_expression[sel_genes, ]
    }
    
    ## compute the correlation of TF activity and gene expression along the trajectory
    ## df.cor -> gene by TF matrix
    df.cor <- t(cor(t(tf_activity), t(gene_expression)))
    
    if(!is.null(tf.use)){
        df.cor <- df.cor[, tf.use]
    }

    return(df.cor)

}

GetGRN <- function(object, df.cor = NULL, df.p2g = NULL, 
                   tf_assay = "chromvar", gene_assay = "RNA",
                   trajectory.name = "Trajectory", min.cor = 0.4){

    if(is.null(df.cor)){
        stop("Please provide TF activity and expression correlation results!")
    }
    
    if(is.null(df.p2g)){
        stop("Please provide peak-to-gene links!")
    }
    
    ## get tf activity and gene expression along trajectory
    trajMM <- GetTrajectory(object, 
                        assay = tf_assay, 
                        slot = "data",
                            trajectory.name = trajectory.name,
                       smoothWindow = 7,
                       log2Norm = FALSE)

    trajRNA <- GetTrajectory(object, 
                        assay = gene_assay, 
                        slot = "data",
                             trajectory.name = trajectory.name,
                       smoothWindow = 7,
                       log2Norm = TRUE)
    
    rownames(trajMM) <- object@assays$ATAC@motifs@motif.names

    tf_activity <- suppressMessages(TrajectoryHeatmap(trajMM, 
                            varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                           name = "TF activity",
                           returnMatrix = TRUE)) 

    gene_expression <- suppressMessages(TrajectoryHeatmap(trajRNA, 
                            varCutOff = 0.9,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                           name = "Gene expression",
                                        returnMatrix = TRUE))
    
    ## here we filter the TFs according to our correlation analysis
    tf_activity <- tf_activity[df.cor$tfs, ]
    
    ## here we further filter the genes by only consider genes that are linked to peaks
    sel_genes <- intersect(rownames(gene_expression), unique(df.p2g$gene))
    
    ## subset the gene expression matrix
    gene_expression <- gene_expression[sel_genes, ]
    
    ## compute the correlation of TF activity and gene expression along the trajectory
    ## df.grn -> gene by TF matrix
    df.grn <- t(cor(t(tf_activity), t(gene_expression)))
    df.grn <- df.grn[, df.cor$tfs]
    
    ## We next use peak-to-gene links to predict the target genes for each TF.
    ## We consider a gene is regulated by a peak if there is a positive
    ## correlation between gene expression and peak accessibility
    ## mat.p2g is a gene by peak data frame
    mat.p2g <- df.p2g %>%
        select(c(peak, gene, Correlation)) %>%
        tidyr::pivot_wider(names_from = peak, values_from = Correlation) %>%
        textshape::column_to_rownames("gene")

    mat.p2g[is.na(mat.p2g)] <- 0
    mat.p2g[mat.p2g>0] <- 1
    
    ## To link gene to TF, we also need the TF binding information
    ## Here we obtain a peak by TF matrix representing if peak is bound by a TF 
    ## mat.motif is a peak by TF matrix
    mat.motif <- object@assays$ATAC@motifs@data
    colnames(mat.motif) <- object@assays$ATAC@motifs@motif.names
    
    ## We can filter this complete matching matrix by only using the peaks that are
    ## linked with some genes, and TFs that are selected based on correlation analysis
    mat.motif <- mat.motif[colnames(mat.p2g), colnames(df.grn)]
    
    ## multiply gene by peak and peak by tf matrix to obtain gene by tf matrix
    mat.tf.gene <- as.matrix(mat.p2g) %*% as.matrix(mat.motif)
    
    ## we binarize this matrix to indicate if a gene is regulated by a TF through a peak
    mat.tf.gene[mat.tf.gene > 0] <- 1
    
    ## here we multiply correlation matrix and regulation matrix
    df.grn <- as.matrix(df.grn) * mat.tf.gene %>%
        as.data.frame()
    df.grn$gene <- rownames(df.grn)
    
    df.grn <- df.grn %>%
       tidyr::pivot_longer(!gene, names_to = "tf", values_to = "correlation") %>%
       subset(correlation > min.cor)

    return(df.grn)

}

GRNPlot <- function(df.grn){
    g <- igraph::graph_from_data_frame(df.grn, directed = FALSE)
    layout <- igraph::layout_with_fr(g, weights = E(g)$correlation)
    
    ## compute pagerank and betweenness
    pagerank <- page_rank(netobj,weights = E(netobj)$correlation)
    bet <-betweenness(netobj,weights = E(netobj)$correlation,normalized = TRUE)

}