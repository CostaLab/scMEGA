#' Add trajectory
#'
#' This function infers a pseudo-time trajectory based on the supervised input.
#' It was modified from the package ArchR \code{\link{addTrajectory}}
#' to use a Seurat object as input.
#' For more details, check here \url{https://www.archrproject.com/reference/addTrajectory.html}.
#' @param object A Seurat object
#' @param name A string indicating the name of the fitted trajectory to be added in meta.data
#' @param trajectory The order of cell groups to be used for constraining
#' the initial supervised fitting procedure
#' @param group.by The column name in metadata that contains the cell group
#' definitions used in \code{trajectory}
#' @param reduction Name of dimension reduction used to infer the trajectory
#' @param dims How many dimensions to use
#' @param pre.filter.quantile Prior to the initial supervised trajectory fitting,
#' cells whose euclidean distance from the cell-grouping center is
#' above the provided quantile will be excluded.
#' @param post.filter.quantile After initial supervised trajectory fitting, cells
#' whose euclidean distance from the cell-grouping center is above the provided
#' quantile will be excluded.
#' @param use.all Whether or not to use cells ourside of trajectory groups for
#' post-fitting procedure
#' @param dof The number of degrees of freedom to be used in the spline fit
#' @param spar The sparsity to be used in the spline fit.
#' @param seed Random seed
#'
#' @return Seurat object with inferred trajectory
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- AddTrajectory(
#'    object = obj,
#'    trajectory = c(2, 1),
#'    group.by = "RNA_snn_res.0.2",
#'    reduction = "dm",
#'    dims = 1:3,
#'    use.all = FALSE
#' )
#' }
AddTrajectory <-
  function(object,
           name = "Trajectory",
           trajectory = NULL,
           group.by = NULL,
           reduction = NULL,
           dims = NULL,
           pre.filter.quantile = 0.9,
           post.filter.quantile = 0.9,
           use.all = FALSE,
           dof = 250,
           spar = 1,
           seed = 42) {
    set.seed(seed)

    if (is.null(group.by)) {
      stop("Please specify the group based on which the trajectory is inferred!")
    }

    if (is.null(trajectory)) {
      stop("Please specify the trajectory!")
    }

    if (is.null(reduction)) {
      stop("Please prodive the dimensional reduction!")
    }

    #    df.group <- DataFrame(object@meta.data[, group.by])
    df.group <- object@meta.data[, group.by] %>%
      as.data.frame()
    rownames(df.group) <- colnames(object)
    df.group <- df.group[df.group[, 1] %in% trajectory, , drop = F]

    if (sum(unique(df.group[, 1]) %in% trajectory) == 0) {
      stop("Cannot find the specified trajectory!")
    }

    if (is.null(dims)) {
      data.use <- Seurat::Embeddings(object, reduction = reduction)
    } else{
      data.use <-
        Seurat::Embeddings(object, reduction = reduction)[, dims]
    }

    data.use <- data.use[rownames(df.group), , drop = FALSE]

    ######################################################
    #Filter Outliers
    ######################################################

    filterObj <- lapply(seq_along(trajectory), function(x) {
      #Subset
      groupsx <- rownames(df.group)[df.group[, 1] == trajectory[x]]
      matx <- data.use[groupsx, , drop = FALSE]

      #Filter Distance
      matMeanx <- colMeans(matx)
      diffx <- sqrt(colSums((t(matx) - matMeanx) ^ 2))
      idxKeep <-
        which(diffx <= quantile(diffx, pre.filter.quantile))

      #Filter
      list(mat = matx[idxKeep, , drop = FALSE],
           groups = groupsx[idxKeep])

    })

    matFilter <-
      lapply(seq_along(filterObj), function(x)
        filterObj[[x]]$mat) %>% Reduce("rbind", .)
    groupsFilter <-
      df.group[lapply(seq_along(filterObj), function(x)
        filterObj[[x]]$groups) %>% Reduce("c", .), , drop = FALSE]

    ######################################################
    #Now Initial Alignment
    ######################################################

    initialTime <- lapply(seq_along(trajectory), function(x) {
      groupsx <-
        rownames(groupsFilter)[groupsFilter[, 1] == trajectory[x]]
      matx <- matFilter[groupsx, , drop = FALSE]

      #Get Differences
      if (x != length(trajectory)) {
        groupsxp1 <-
          rownames(groupsFilter)[groupsFilter[, 1] == trajectory[x + 1]]
        meanx <- colMeans(matFilter[groupsxp1, , drop = FALSE])
        diffx <- sqrt(colSums((t(matx) - meanx) ^ 2))
        timex <- (1 - ArchR:::.getQuantiles(diffx)) + x
      } else{
        groupsxm1 <-
          rownames(groupsFilter)[groupsFilter[, 1] == trajectory[x - 1]]
        meanx <- colMeans(matFilter[groupsxm1, , drop = FALSE])
        diffx <- sqrt(colSums((t(matx) - meanx) ^ 2))
        timex <- ArchR:::.getQuantiles(diffx) + x
      }

      timex

    }) %>% unlist

    ######################################################
    #Fit Cubic Splines
    ######################################################

    matSpline <- lapply(seq_len(ncol(matFilter)), function(x) {
      tryCatch({
        stats::smooth.spline(
          x = initialTime,
          y = matFilter[names(initialTime), x],
          df = dof,
          spar = spar
        )[[2]]
      }, error = function(e) {
        errorList <- list(
          it = x,
          x = initialTime,
          y = matFilter[names(initialTime), x],
          df = dof,
          spar = spar
        )

      })
    }) %>% Reduce("cbind", .) %>% data.frame()

    ######################################################
    # 1. KNN Fit vs Actual
    ######################################################

    knnObj <- nabor::knn(data =  matSpline,
                         query = data.use,
                         k = 3)

    #Estimate place along trajectory
    knnIdx <- knnObj[[1]]
    knnDist <- knnObj[[2]]
    knnDiff <- ifelse(knnIdx[, 2] > knnIdx[, 3], 1,-1)
    knnDistQ <- ArchR:::.getQuantiles(knnDist[, 1])

    #Filter Outlier Cells to Trajectory for High Resolution
    idxKeep <-
      which(knnDist[, 1] <= quantile(knnDist[, 1], post.filter.quantile))
    dfTrajectory <- data.frame(
      row.names = rownames(data.use),
      Distance = knnDist[, 1],
      DistanceIdx = knnIdx[, 1] + knnDiff * knnDistQ
    )[idxKeep, , drop = FALSE]


    dfTrajectory3 <- dfTrajectory

    dfTrajectory3$Trajectory <-
      100 * ArchR:::.getQuantiles(dfTrajectory3[, 2])
    nas <- rep(NA, dim(object)[2])
    names(nas) <- colnames(object)
    nas[rownames(dfTrajectory3)] <- dfTrajectory3$Trajectory
    object@meta.data[, name] <- nas

    return(object)


  }


#' Get data along the trajectory
#'
#' This function collects and smooths data along the trajectory using Seurat
#' object as input.
#' It was modified from the package ArchR \code{\link{getTrajectory}}
#' to use a Seurat object as input.
#' For more details, check here \url{https://www.archrproject.com/reference/getTrajectory.html}.
#' @param object Seurat object
#' @param trajectory.name The name of trajectory inferred by \code{\link{AddTrajectory}}
#' @param assay Which assay is used for data collection
#' @param slot Name of slot used to collect data. Usually set as "counts" or "data".
#' @param groupEvery The number of sequential percentiles to group together when generating a trajectory.
#' This is similar to smoothing via a non-overlapping sliding window across pseudo-time.
#' @param log2Norm Whether or not the data should be log2 transformed. This can be set to TRUE
#' if gene expression or chromatin accessibility data is used,
#' otherwise FALSE if data from TF activity is used.
#' @param scaleTo Once the sequential trajectory matrix is created,
#' each column in that matrix will be normalized to a column sum indicated by scaleTo.
#' @param smoothWindow An integer value indicating the smoothing window in size
#' (relative to groupEvery) for the sequential trajectory matrix to better reveal
#' temporal dynamics.
#'
#' @import Seurat
#' @return A SummarizedExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' trajRNA <- GetTrajectory(
#'    object = obj,
#'    assay = "RNA",
#'    slot = "data",
#'    smoothWindow = 7,
#'    log2Norm = TRUE
#' )
#' }
GetTrajectory <- function(object = NULL,
                          trajectory.name = "Trajectory",
                          assay = NULL,
                          slot = "counts",
                          groupEvery = 1,
                          log2Norm = TRUE,
                          scaleTo = 10000,
                          smoothWindow = 11) {
  if (is.null(assay) | !assay %in% Assays(object)) {
    stop("Please provide an available assay!")
  }

  if (!(trajectory.name %in% colnames(object@meta.data))) {
    stop(glue::glue("Cannot find trajecotry {trajectory.name}!"))
  }

  trajectory <- object@meta.data[trajectory.name]

  trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]
  breaks <- seq(0, 100, groupEvery)
  if (!all(is.numeric(trajectory[, 1]))) {
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if (!all(trajectory[, 1] >= 0 & trajectory[, 1] <= 100)) {
    stop(
      "Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?"
    )
  }

  groupList <- lapply(seq_along(breaks), function(x) {
    if (x == 1) {
      NULL
    } else{
      rownames(trajectory)[which(trajectory[, 1] > breaks[x - 1] &
                                   trajectory[, 1] <= breaks[x])]
    }
  })[-1]

  names(groupList) <-
    paste0("T.", breaks[-length(breaks)], "_", breaks[-1])

  message("Creating Trajectory Group Matrix..")
  data.use <- GetAssayData(object, assay = assay, slot = slot)

  groupMat <- lapply(1:length(groupList), function(x) {
    cell_names <- groupList[[x]]
    mat <- Matrix::rowMeans(data.use[, cell_names])

  }) %>% Reduce(cbind, .)

  colnames(groupMat) <- names(groupList)

  #Scale
  if (!is.null(scaleTo)) {
    if (any(groupMat < 0)) {
      message(
        "Some values are below 0, this could be the Motif activity matrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!"
      )
    } else{
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }
  }

  if (log2Norm) {
    if (any(groupMat < 0)) {
      message(
        "Some values are below 0, this could be a Motif activity matrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!"
      )
    } else{
      groupMat <- log2(groupMat + 1)
    }
  }

  if (!is.null(smoothWindow)) {
    message("Smoothing...")
    smoothGroupMat <-
      as.matrix(t(apply(groupMat, 1, function(x)
        ArchR:::.centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(assays = SimpleList(
      smoothMat = as.matrix(smoothGroupMat),
      mat = as.matrix(groupMat)
    ))

  } else{
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <-
      SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)))

  }

  return(seTrajectory)
}


#' Get correlation
#'
#' This function will compute the correlation between TF expression as measured
#' by RNA-seq and binding activity estimated from ATAC-seq along the trajectory.
#'
#' @param trajectory1 A SummarizedExperiment object for TF activity
#' generated by the function \code{\link{GetTrajectory}}
#' @param trajectory2 A SummarizedExperiment object for TF expression
#' generated by the function \code{\link{GetTrajectory}}
#'
#' @return A data frame cantaining the correlation
#' @export
#'
#' @examples
#' \dontrun{
#' df.cor <- GetCorrelation(
#'    trajectory1 = trajMM,
#'    trajectory2 = trajRNA
#' )
#' }
GetCorrelation <- function(trajectory1,
                           trajectory2) {
  mat1 <- assay(trajectory1)
  mat2 <- assay(trajectory2)

  features.use <- intersect(rownames(mat1), rownames(mat2))

  if (length(features.use) == 0) {
    stop("Cannot find any shared features!")
  } else{
    message(glue::glue("Find {length(features.use)} shared features!"))
  }

  mat1 <- mat1[features.use,]
  mat2 <- mat2[features.use,]

  df_cor <- lapply(1:length(features.use), function(x) {
    cor.res <- suppressWarnings(cor.test(mat1[x,], mat2[x,]))
    p.value <- cor.res$p.value
    estimate <- cor.res$estimate[[1]]
    df <- data.frame(estimate, p.value)
    rownames(df) <- features.use[[x]]

    colnames(df) <- c("correlation", "p_value")

    return(df)

  }) %>% Reduce(rbind, .) %>%
    tidyr::drop_na()

  df_cor$adj_p <- p.adjust(df_cor$p_value, method = "fdr")
  df_cor$tfs <- rownames(df_cor)

  return(df_cor)

}
