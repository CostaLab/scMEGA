GetTFGeneCorrelation <- function(object, 
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

GetGRN <- function(object,
                   df.cor = NULL, 
                   df.p2g = NULL,
                   min.cor = 0.4){

    if(is.null(df.cor)){
        stop("Please provide a tf-gene correlation matrix!")
    }
    
    if(is.null(df.p2g)){
        stop("Please provide peak-to-gene links!")
    }
    
    ## We next use peak-to-gene links to predict the target genes for each TF.
    ## We consider a gene is regulated by a peak if there is a positive
    ## correlation between gene expression and peak accessibility
    ## mat.p2g is a gene by peak data frame
    message("Filtering network by peak-to-gene links...")
    mat.p2g <- df.p2g %>%
        select(c(peak, gene, Correlation)) %>%
        tidyr::pivot_wider(names_from = peak, values_from = Correlation) %>%
        textshape::column_to_rownames("gene")

    mat.p2g[is.na(mat.p2g)] <- 0
    mat.p2g[mat.p2g>0] <- 1
    
    ## To link gene to TF, we also need the TF binding information
    ## Here we obtain a peak by TF matrix representing if peak is bound by a TF 
    ## mat.motif is a peak by TF matrix
    message("Filtering network by TF binding site...")
    mat.motif <- object@assays$ATAC@motifs@data
    colnames(mat.motif) <- object@assays$ATAC@motifs@motif.names
    
    ## We can filter this complete matching matrix by only using the peaks that are
    ## linked with some genes, and TFs that are selected based on correlation analysis
    mat.motif <- mat.motif[colnames(mat.p2g), colnames(df.cor)]
    
    ## multiply gene by peak and peak by tf matrix to obtain gene by tf matrix
    mat.tf.gene <- as.matrix(mat.p2g) %*% as.matrix(mat.motif)
    
    ## we binarize this matrix to indicate if a gene is regulated by a TF through a peak
    mat.tf.gene[mat.tf.gene > 0] <- 1
    
    ## here we multiply correlation matrix and regulation matrix
    df.grn <- as.matrix(df.cor) * mat.tf.gene %>%
        as.data.frame()
    
    df.grn$gene <- rownames(df.grn)
    
    df.grn <- df.grn %>%
       tidyr::pivot_longer(!gene, names_to = "tf", values_to = "correlation") %>%
       subset(correlation > min.cor)

    return(df.grn)

}

GRNHeatmap <- function(tf.gene.cor, tf.timepoint, km = 2){
    col_fun <- circlize::colorRamp2(tf.timepoint, 
                                    ArchR::paletteContinuous(set = "blueYellow", 
                                                             n = length(tf.timepoint)))

    column_ha <- ComplexHeatmap::HeatmapAnnotation(time_point = tf.timepoint,
                                                   col = list(time_point = col_fun))

    ht <- Heatmap(as.matrix(tf.gene.cor),
                   name = "correlation",
                   cluster_columns = FALSE,
                 clustering_method_rows = "ward.D2",
                   top_annotation = column_ha,
                  show_row_names = FALSE,
                   show_column_names = TRUE,
                    row_km = km,
                  column_km = km,
                  border = TRUE)
    
    return(ht)


}

GRNPlot <- function(df.grn,
                    tfs.timepoint = NULL, 
                    genes.cluster = NULL,
                    genes.highlight = NULL,
                    cols.highlight = "#984ea3", seed = 42){
    
    if(is.null(tfs.timepoint)){
        stop("Need time point for each TF!")
    }
    
    # create graph from data frame
    g <- graph_from_data_frame(df.grn, directed=TRUE)
   
    # compute pagerank and betweenness
    pagerank <- page_rank(g, weights = E(g)$correlation)
    bet <-betweenness(g, weights = E(g)$correlation, normalized = TRUE)
    df_measure <- data.frame(node = V(g)$name, 
                              pagerank = pagerank$vector, 
                              betweenness = bet) %>%
        subset(node %in% df.grn$tf)
    
    # compute importance only for TFs based on centrality and betweenness
    min.page <- min(df_measure$pagerank)
    min.bet <- min(df_measure$betweenness)
    df_measure$importance <- sqrt((df_measure$pagerank - min.page)**2 + 
                                 (df_measure$betweenness - min.bet)**2)
    
    # assign size to each node
    # for TFs, the size is proportional to the importance
    tf_size <- df_measure$importance
    
    ## for genes, we use the minimum size of TFs
    gene_size <- rep(min(df_measure$importance), length(unique(df.grn$gene)))
    
    # assign color to each node
    ## TFs are colored by pseudotime point
    cols.tf <- ArchR::paletteContinuous(set = "blueYellow", 
                                         n = length(tfs.timepoint))
    names(cols.tf) <- names(tfs.timepoint)
    
    ## genes are colored based on the clustering
    cols.gene <- ArchR::paletteDiscrete(values = names(genes.cluster))

    names(gene_color) <- df_gene_clustering$gene

    v_color <- c(tf_color, gene_color)
    
    
    # compute layout
    set.seed(seed)
    layout <- layout_with_fr(g, weights = E(g)$correlation, dim = 2, niter = 1000)
    
    p <- ggraph(g, layout = layout) +
    geom_edge_link(edge_colour = "gray", edge_alpha = 0.25)
    
    
    # highlight some TFs
    if(!is.null(genes.highlight)){
        p <- p + geom_node_label(aes(filter = V(g)$name %in% genes.highlight,
                    label = V(g)$name),
                    repel = TRUE,
                    hjust = "inward",
                size = 5,
                color = cols.highlight,
                show.legend = FALSE)
    
    }

    p <- p + theme_void()
    
    return(p)
}