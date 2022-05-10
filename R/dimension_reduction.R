#' Run diffusion map
#'
#' This function is a wrapper to run the diffusion map dimensional reduction
#' technique via package destiny. Details on this package can be found here:
#'\url{https://bioconductor.org/packages/release/bioc/html/destiny.html}.
#' @param object A Seurat object
#' @param reduction.key dimensional reduction key,
#' specifies the string before the number for the dimension names. DC by default
#' @param assay Assay to pull data
#' @param reduction Which dimensional reduction to use for the diffusion map input
#' @param dims Which dimensions to use as input features
#' @param k Number of nearest neighbors to consider. This is passed to the function
#' DiffusionMap
#' @param seed Random seed
#' @param verbose Whether or not print running message
#' @param reduction.name Name to store dimensional reduction under in the Seurat object.
#' Default: dm
#'
#' @import destiny
#' @import Seurat
#' @return Returns a Seurat object containing a diffusion map representation
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' # Run diffusion map on first 5 PCs
#' pbmc_small <- RunDiffusionMap(object = pbmc_small, reduction = "pca",
#' dims = 1:20, k = 10)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction = 'dm')
#' }
RunDiffusionMap <-
  function(object,
           reduction.key = "DC_",
           assay = DefaultAssay(object = object),
           reduction = NULL,
           dims = NULL,
           k = NULL,
           seed = 42,
           verbose = TRUE,
           reduction.name = "dm") {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    if (is.null(reduction)) {
      stop("Please specify the dimension reduction based on which diffusion map is performed!")
    }


    ## diffusion map
    if (!is.null(dims)) {
      data.use <- Embeddings(object[[reduction]])[, dims]
    } else{
      data.use <- Embeddings(object[[reduction]])
    }

    if (!is.null(k)) {
      dm <- DiffusionMap(as.matrix(data.use),
                         verbose = verbose,
                         k = 30)
    } else{
      dm <- DiffusionMap(as.matrix(data.use),
                         verbose = verbose)
    }

    dm <- as.data.frame(dm)
    colnames(dm) <- paste0(reduction.key, 1:ncol(dm))

    object[[reduction.name]] <-
      CreateDimReducObject(embeddings = as.matrix(dm),
                           key = reduction.key,
                           assay = assay)

    return(object)


  }
