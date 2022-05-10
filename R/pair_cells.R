CoembeddingData <- function(obj.rna, obj.atac, gene.activity, reference.assay = "RNA",
                       reduction = "cca", weight.reduction = NULL, verbose = TRUE){
    ## make sure there are the same number of cells in atac and gene activity data
    if(ncol(obj.atac) != ncol(gene.activity)){
        stop("The number of cells in ATAC-seq and Gene activity are not the same!")
    }
    
    if(is.null(weight.reduction)){
        stop("Please specify dimensional reduction to use for the weighting anchors!")
    }
    
    gene.use <- intersect(rownames(gene.activity),
                          rownames(obj.rna))

    message("Performing data integration using Seurat...")
    
    obj.atac[['RNA']] <- CreateAssayObject(counts = gene.activity[gene.use, ])
    DefaultAssay(obj.atac) <- "RNA"
    
    obj.atac <- obj.atac %>%
        NormalizeData(verbose = verbose) %>%
        FindVariableFeatures(verbose = verbose) %>%
        ScaleData(verbose = verbose)
    
    transfer.anchors <- FindTransferAnchors(reference = obj.rna,
                                          query = obj.atac,
                                          features = gene.use,
                                          reference.assay = reference.assay,
                                          query.assay = "RNA",
                                          reduction = reduction,
                                           verbose = verbose)
    
    # we here restrict the imputation to the selected genes
    refdata <- GetAssayData(obj.rna, assay = reference.assay, slot = "data")[gene.use, ]

    # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells. 
    # imputation (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors,
                               refdata = refdata,
                               weight.reduction = obj.atac[[weight.reduction]],
                               dims = 1:30)

    obj.imputation <- CreateSeuratObject(counts = imputation, assay = "RNA")
    DefaultAssay(obj.imputation) <- "RNA"
    
    meta.data.use <- intersect(colnames(obj.atac@meta.data),
                              colnames(obj.rna@meta.data))
    
    obj.imputation <- AddMetaData(obj.imputation, obj.atac@meta.data[, meta.data.use])
    obj.rna@meta.data <- obj.rna@meta.data[, meta.data.use]
    
    obj.rna$tech <- "RNA"
    obj.imputation$tech <- "ATAC"
    
    # merge the objects
    coembed <- merge(x = obj.imputation, y = obj.rna)

    #Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both datasets
    message("Coemebdding the data...")
    coembed <- coembed %>%
        ScaleData(features = gene.use, do.scale = FALSE) %>%
        FindVariableFeatures() %>%
        RunPCA(verbose = FALSE) %>%
        RunUMAP(dims = 1:30, verbose = FALSE)
    
    return(coembed)
}


PairingCells <- function(obj.coembed, reduction = NULL, graph.name = NULL,
                         paired.by = NULL,
                         assay = "RNA", mode = "geodesic",
                         tol = 0.0001, search.range = 0.2, max.multimatch = 5, 
                         min_subgraph_size = 50,seed = 42,
                        k = 300, n.cores = 30){
    
    if(is.null(reduction)){
        stop("Please specify dimensional reduction to use for the pairing cells!")
    }
    
    message("Getting dimensional reduction data for pairing cells...")
    obj.atac <- subset(obj.coembed, tech == "ATAC")
    obj.rna <- subset(obj.coembed, tech == "RNA")
    
    embedding.atac <- Embeddings(object = obj.atac, reduction = reduction)
    embedding.rna <- Embeddings(object = obj.rna, reduction = reduction)
    
    embedding <- rbind(embedding.atac, embedding.rna)
    n.cells <- dim(embedding)[1]
    
    # release memory
    rm(obj.atac)
    rm(obj.rna)
    gc()
    
    message(glue("Pairing cells using {mode} mode..."))
    
    if (mode == "geodesic"){
        message("Constructing KNN graph for computing geodesic distance ..")
        
        obj.coembed <- FindNeighbors(obj.coembed, reduction = reduction, assay = "RNA",
                                    verbose = FALSE)
        
        adjmatrix <- obj.coembed@graphs$RNA_nn
        diag(adjmatrix) <- 0
        #adj.matrix <- as.matrix(x = obj.coembed@graphs$RNA_nn)
        #adj.matrix <- adj.matrix + t(adj.matrix)
        knn.graph <- graph_from_adjacency_matrix(adjmatrix = adjmatrix,
                                                 mode = "max",
                                                weighted = NULL)
        
        message("Computing graph-based geodesic distance ..")
        # Compute shortest paths between all cell pairs. 
        # We checked that this returns a symmetric matrix
        shortest.paths <- shortest.paths(knn.graph)

        # Find connected subgraphs
        sub.graphs <- clusters(knn.graph)
        message(glue("# KNN subgraphs detected: {length(unique(sub.graphs$membership))}"))
        
        obj.coembed$subgraph <- sub.graphs$membership

        all.pairs <- NULL

        message("Skipping subgraphs with either ATAC/RNA cells fewer than: ", 
                min_subgraph_size)

        # Go through each subgraph
        for (sub.graph.idx in unique(sub.graphs$membership)){
          message("Pairing cells for subgraph No.", sub.graph.idx)

          # Retrieve the subgraph
          sub.graph.nodes <- sub.graphs$membership == sub.graph.idx
          knn.sub.graph <- induced_subgraph(knn.graph, which(sub.graph.nodes))

          # Use down-sampling to make sure in this subgraph the number of ATAC and RNA cells are balanced
          subgraph_cells <- colnames(obj.coembed)[sub.graph.nodes]
          n_ATAC <- sum(subgraph_cells %in% rownames(embedding.atac))
          n_RNA <- sum(subgraph_cells %in% rownames(embedding.rna))

          message("Total ATAC cells in subgraph: ",n_ATAC)
          message("Total RNA cells in subgraph: ",n_RNA)

          embedding.atac.sub <- embedding.atac[sub.graph.nodes[1: dim(embedding.atac)[1]],]
          embedding.rna.sub <- embedding.rna[sub.graph.nodes[(dim(embedding.atac)[1] + 1):n.cells], ]

          if (n_ATAC > n_RNA){
            set.seed(seed)
            embedding.atac.sub <- embedding.atac.sub[sample(1:n_ATAC, n_RNA, replace = FALSE),]
          } else if (n_ATAC < n_RNA){
            set.seed(seed)
            embedding.rna.sub <- embedding.rna.sub[sample(1:n_RNA, n_ATAC, replace = FALSE),]
          }

          if(is.null(nrow(embedding.atac.sub)) | is.null(nrow(embedding.rna.sub))){
            message("Down-sampling within subgraph between assays led to 0 cells in one assay..")
            message("Skipping current subgraph")
            next
          }

      # Subset the original geodesic distance matrix to get the geodesic distance matrix for the subgraph
      subgraph_geodist <- shortest.paths[match(rownames(embedding.atac.sub), rownames(embedding)),
                                         match(rownames(embedding.rna.sub), rownames(embedding))]
      subgraph_size <- dim(subgraph_geodist)[1]

      message("Subgraph size: ", subgraph_size)

      # TO AVOID MAJOR SUBGRAPH(S) WERE BEING SKIPPED SOMETIMES
      size_threshold <- ceiling((nrow(embedding.atac.sub) + nrow(embedding.rna.sub)) * search.range)
      k_pairing <- size_threshold

      message("Search threshold being used: ",k_pairing)

      if (subgraph_size < size_threshold | subgraph_size < min_subgraph_size) {
        message("Insufficient number of cells in subgraph. Skipping current subgraph")
        next
      }

      # We also calculate euclidean distance matrix
      subgraph_eucdist <- distmat(embedding.atac.sub, embedding.rna.sub)

      # Find KNN based on geodesic distances.
      message("Constructing KNN based on geodesic distance to reduce search pairing search space")
      geodist_knn <- array(-1, dim = dim(subgraph_geodist))
            
      for (i in 1:subgraph_size){
        # Find RNA cells in the KNN of each ATAC cell
        geodist_threshold <- sort(subgraph_geodist[i, ])[k_pairing]
        knn_ind <- subgraph_geodist[i,] < geodist_threshold
        geodist_knn[i, knn_ind] <- subgraph_eucdist[i, knn_ind]

        # Find ATAC cells in the KNN of each RNA cell
        geodist_threshold <- sort(subgraph_geodist[, i])[k_pairing]
        knn_ind <- subgraph_geodist[,i] < geodist_threshold
        geodist_knn[knn_ind, i] <- subgraph_eucdist[knn_ind, i]

      }

      # For an ATAC-RNA cell pair, if neither of them are in the KNN of the other, we set their distance to be inf
      geodist_knn[geodist_knn < 0] <- Inf

      # Add rownames to the matrix
      rownames(geodist_knn) <- paste0("ATAC_", 1:subgraph_size)
      colnames(geodist_knn) <- paste0("RNA_", 1:subgraph_size)

      message(paste("Number of cells being paired:",
                  dim(geodist_knn)[1],
                  "ATAC and",
                  dim(geodist_knn)[1],
                  " RNA cells"))

      message("Determing pairs through optimized bipartite matching ..\n")
      options("optmatch_max_problem_size" = Inf)
            
      cell_matches <- suppressWarnings(fullmatch(as.InfinitySparseMatrix(as.matrix(geodist_knn)),
                                                           tol = tol,
                                                           min.controls = 1 / max.multimatch,
                                                           max.controls = max.multimatch))
            
      pair.list <- GetPairList(cell_matches,
                                 rownames(embedding.atac.sub),
                                 rownames(embedding.rna.sub))

      message("Finished!\n")
      
      # Append the results for this subgraph to the list of all results
      all.pairs <- rbind(all.pairs, pair.list)
    }
  }
    
    # pairs are sometimes repreated, here we make the results unique
#    all.pairs <- unique(all.pairs)
    
    #print(length(unique(all.pairs)))
    all.pairs <- all.pairs[!duplicated(all.pairs$ATAC), ]
    all.pairs$cell_name <- paste0("cell_", 1:nrow(all.pairs))
    
    return(all.pairs)
}

CreatePairedObject <- function(df.pair, obj.rna, obj.atac, use.assay1 = NULL, use.assay2 = NULL,
                              sep = c("-", "-")){
    
    if(is.null(use.assay1)){
        stop("Please provide the name for assay from RNA object")
    }
    
    if(is.null(use.assay2)){
        stop("Please provide the name for assay from ATAC object")
    }
    
    message("Merging the count data...")
    rna.counts <- GetAssayData(obj.rna, assay = use.assay1, slot = "counts")[, df.pair$RNA]
    atac.counts <- GetAssayData(obj.atac, assay = use.assay2, slot = "counts")[, df.pair$ATAC]

    rna.data <- GetAssayData(obj.rna, assay = use.assay1, slot = "data")[, df.pair$RNA]
    atac.data <- GetAssayData(obj.atac, assay = use.assay2, slot = "data")[, df.pair$ATAC]
    
    colnames(rna.counts) <- df.pair$cell_name
    colnames(rna.data) <- df.pair$cell_name
    
    colnames(atac.counts) <- df.pair$cell_name
    colnames(atac.data) <- df.pair$cell_name
    
    # create a Seurat object containing the RNA adata
    obj <- CreateSeuratObject(
      counts = rna.counts,
        data = rna.data,
      assay = use.assay1
    )

    # create ATAC assay and add it to the object
    obj[[use.assay2]] <- CreateChromatinAssay(
      counts = atac.counts,
      sep = sep,
        min.cells = 10
    )

    message("Merging the dimension reductions...")
    DefaultAssay(obj) <- use.assay1
    for(reduction in names(obj.rna@reductions)){
        embedding <- Embeddings(obj.rna, reduction = reduction)[df.pair$RNA, ]
        rownames(embedding) <- df.pair$cell_name
        obj[[glue("{reduction}_{use.assay1}")]] <- CreateDimReducObject(embeddings = embedding,  
                                                                       assay = DefaultAssay(obj))
    }
    
    DefaultAssay(obj) <- use.assay2
    for(reduction in names(obj.atac@reductions)){
        embedding <- Embeddings(obj.atac, reduction = reduction)[df.pair$ATAC, ]
        rownames(embedding) <- df.pair$cell_name
        obj[[glue("{reduction}_{use.assay2}")]] <- CreateDimReducObject(embeddings = embedding,  
                                                                       assay = DefaultAssay(obj))
    }
      
    # add metadata, here we use the metadata from RNA assay
    meta.data <- obj.rna@meta.data[df.pair$RNA, ]
    rownames(meta.data) <- df.pair$cell_name

    obj <- AddMetaData(obj, metadata = meta.data)
    
    return(obj)

}


# This function takes in the output of the fullmatch function and sort the results into a list of ATAC-RNA pairs
GetPairList <- function(cell_matches, # The output object of fullmatch or pairmatch in the package optmatch
                          ATAC_barcodes, # Barcode of ATAC cells, must match the order of IDs in cell_matches.
                          #e.g. ATAC_1 in cell matches should correspond to the first element in ATAC_barcodes
                          RNA_barcodes # Barcode of RNA cells, must match the order of IDs in cell_matches.
                          #e.g. RNA_1 in cell matches should correspond to the first element in RNA_barcodes
){

  cell_matches <- sort(cell_matches) # Sort to get ATAC, RNA tuples

  if(length(cell_matches)==0){
      stop("Matches could not be found .. Perhaps try adjusting the constraints to allow optimal matching to be solved?\n")
  }


  if(any(is.na(cell_matches))){
      warning("NA pairs exist ..\n")
  }

  # Currently the result list contain cells paired to multiple other cells
  # If a cell is paired to k cells, we duplicate this cell k times to make k pairs
  # Thus we generate a new list consisting pairs of 1 ATAC - 1 RNA cells
  # We also make sure that in a pair, the first ID is ATAC and the second ID is RNA
  matches <- character()
  pair_ids <- unique(unname(cell_matches))
  for (pair_id in 1:length(pair_ids)){
    new_match <- names(cell_matches[unname(cell_matches)== pair_ids[pair_id]])
    new_match_ATAC <- new_match[splitAndFetch(new_match,"_",1) == "ATAC"]
    new_match_RNA <- new_match[splitAndFetch(new_match,"_",1) == "RNA" ]
    new_match <- vector()
    if (length(new_match_ATAC) > length(new_match_RNA)){
      new_match[seq(1, 2*length(new_match_ATAC), 2)] <- new_match_ATAC
      new_match[seq(2, 2*length(new_match_ATAC), 2)] <- rep(new_match_RNA, length(new_match_ATAC))
    } else {
      new_match[seq(1, 2*length(new_match_RNA), 2)] <- rep(new_match_ATAC, length(new_match_RNA))
      new_match[seq(2, 2*length(new_match_RNA), 2)] <- new_match_RNA
    }
    matches <- c(matches, new_match)
  }

  # Make sure pair groupings (factors) are adjacent
  #stopifnot(all.equal(cell_matches[seq(1,length(cell_matches),2)],cell_matches[seq(2,length(cell_matches),2)],check.attributes=FALSE))

  ATAC_IDs <- matches[seq(1,length(matches),2)] # Names of ATAC cells
  RNA_IDs <- matches[seq(2,length(matches),2)] # Names of RNA cells

  # Checking if ATAC and RNA tuples are concordant
  stopifnot(all(splitAndFetch(ATAC_IDs,"_",1) %in% "ATAC"))
  stopifnot(all(splitAndFetch(RNA_IDs,"_",1) %in% "RNA"))

  # This is just to make sure 1-1, with the names we gave
  # (can still comprise actual doublets from upsampling if any)
  # stopifnot(all(unique(unique(ATACp))) & all(unique(RNAp)))

  # Get corresponding index relative to input matrix order
  ATAC_inds <- as.numeric(splitAndFetch(ATAC_IDs,"_",2))
  RNA_inds <- as.numeric(splitAndFetch(RNA_IDs,"_",2))

  matches_mat <- matrix(c(ATAC_inds,RNA_inds),ncol=2, byrow = FALSE) # ATAC col1, RNA col2

  message("Assembling pair list ..")
  # Make data frame of matches

  pair_df <- data.frame("ATAC"=ATAC_inds,
                        "RNA"=RNA_inds)

  pair_df <- pair_df %>% arrange(ATAC)

  # Convert to labels if they exist
  pair_df$ATAC <- ATAC_barcodes[pair_df$ATAC] # ATAC cell labels
  pair_df$RNA <- RNA_barcodes[pair_df$RNA] # RNA cell labels

  pair_df
}

splitAndFetch <- function(vec,
                          delim,
                          part){
  if(length(part)==1){
    sapply(strsplit(as.character(vec),delim,fixed=TRUE),"[[",part) } else {
      sapply(strsplit(as.character(vec),delim,fixed = TRUE),function(x) paste(x[part],collapse = delim))
    }
}