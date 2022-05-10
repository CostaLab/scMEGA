RunDiffusionMap <- function(object, reduction.key = "DC_", assay = DefaultAssay(object = object),
                            reduction = NULL, dims = NULL, k = NULL,
                            seed = 42, verbose = TRUE, reduction.name = "dm"
                           ){
    
      if (!is.null(seed)) {
        set.seed(seed)
      }
    
    if(is.null(reduction)){
        stop("Please specify the dimension reduction based on which diffusion map is performed!")
    }
    
    
    ## diffusion map
    if(!is.null(dims)){
        data.use <- Embeddings(object[[reduction]])[, dims]
    } else{
        data.use <- Embeddings(object[[reduction]])
    }
    
    if(!is.null(k)){
        dm <- DiffusionMap(as.matrix(data.use),
                       verbose = verbose,
                      k = 30)
    } else{
        dm <- DiffusionMap(as.matrix(data.use),
                       verbose = verbose)    
    }

    dm <- as.data.frame(dm)
    colnames(dm) <- paste0(reduction.key, 1:ncol(dm))

    object[[reduction.name]] <- CreateDimReducObject(embeddings = as.matrix(dm), 
                                                     key = reduction.key,
                                                    assay = assay)
    
    return(object)


}