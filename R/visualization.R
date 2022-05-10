#' Plot cell proportion
#'
#' This function will generate a bar plot to visualize the cell proportion across
#' different samples
#' @param object Seurat object
#' @param group.by Name of one metadata column to group the cells
#' @param prop.in Name of one metadata column to compute the cell proportion
#' @param cols Specific colors for plotting
#' @import dplyr
#' @import ggplot2
#' @return A ggplot object
#' @export
#'
#' @examples
#' library(Seurat)
#' data(pbmc_small)
#' CellPropPlot(
#' pbmc_small,
#' group.by = "RNA_snn_res.1",
#' prop.in = "groups"
#' )
CellPropPlot <- function(object,
                         group.by = NULL,
                         prop.in = NULL,
                         cols = NULL) {
  if (is.null(group.by)) {
    stop("Please specify how to group the cells!")
  }

  if (is.null(prop.in)) {
    stop("Please specify how to compute cell distribution!")
  }

  df <- object@meta.data %>%
    as.data.frame() %>%
    subset(., select = c(prop.in, group.by)) %>%
    group_by(across(all_of(c(prop.in, group.by)))) %>%
    summarise(counts = n()) %>%
    mutate(proportion = counts / sum(counts))

  proportion <- "proportion"

  p <-
    ggplot(data = df, aes_string(x = prop.in, y = proportion, fill = group.by)) +
    geom_bar(stat = "identity") +
    xlab("") + ylab(proportion) +
    ggtitle(group.by) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  if (!is.null(cols)) {
    p <- p + scale_fill_manual(values = cols)
  }

  return(p)

}

#' Compare the cell proportion
#'
#'
#' @param object Seurat object
#' @param group.by Name of one metadata column to group the cells, usually this refers to
#' the clustering results.
#' @param prop.in Name of one metadata column to compute the cell proportion, usually
#' this refers to the sample identify
#' @param sample.annotation Name of one metadata column to group the samples for comparison
#' @param cols Specific colors for plotting
#' @param comparisons A list of length-2 vectors used to compare the proportion.
#' This parameter is passed to the function stat_compare_means.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- CompareCellProp(
#'    object = coembed,
#'    group.by = "RNA_snn_res.0.2",
#'    prop.in = "patient_region_id",
#'    sample.annotation = "patient_group",
#'     comparisons = list(c("myogenic", "ischemic"),
#'                        c("ischemic", "fibrotic"),
#'                        c("myogenic", "fibrotic"))
#' )
#' }
CompareCellProp <-
  function(object,
           group.by = NULL,
           prop.in = NULL,
           sample.annotation = NULL,
           cols = NULL,
           comparisons = NULL,
           method.test = "wilcox.test") {
    if (is.null(group.by)) {
      stop("Please specify how to group the cells!")
    }
    if (is.null(prop.in)) {
      stop("Please specify how to compute cell proportion!")
    }

    df_prop <- object@meta.data %>%
      as.data.frame() %>%
      subset(., select = c(prop.in, group.by)) %>%
      group_by(across(all_of(c(prop.in, group.by)))) %>%
      summarise(counts = n()) %>%
      mutate(proportion = counts / sum(counts))

    df_anno <- object@meta.data %>%
      as.data.frame() %>%
      subset(., select = c(prop.in, sample.annotation)) %>%
      unique()

    df <- merge.data.frame(df_prop, df_anno)

    proportion <- "proportion"

    if (!is.null(comparisons)) {
      p <-
        ggplot(data = df, aes_string(x = sample.annotation, y = proportion)) +
        geom_boxplot(aes_string(color = sample.annotation)) +
        facet_wrap(as.formula(paste("~", group.by)), nrow = 1) +
        ggpubr::stat_compare_means(comparisons = comparisons,
                           method = method.test) +
        theme_cowplot() +
        xlab("") + ylab("") +
        theme(axis.text.x = element_blank(),
              legend.title = element_blank())

    } else{
      p <-
        ggplot(data = df, aes_string(x = sample.annotation, y = proportion)) +
        geom_boxplot(aes_string(color = sample.annotation)) +
        facet_wrap(as.formula(paste("~", group.by)), nrow = 1) +
        theme_cowplot() +
        xlab("") + ylab("") +
        theme(axis.text.x = element_blank(),
              legend.title = element_blank())
    }

    return(p)

  }

#' Plot trajectory
#'
#' This function generates a scatter plot to visualize the inferred trajectory.
#' It was modified from the package ArchR \code{\link{plotTrajectory}}
#' to use a Seurat object as input.
#' For more details, check here \url{https://www.archrproject.com/reference/plotTrajectory.html}.
#' @param object Seurat object
#' @param trajectory Name of one metadata column to group the cells, usually this refers to
#' the clustering results
#' @param reduction Which dimension reduction is used to visualize the trajectory
#' @param size A number indicating the size of the points to plot
#' @param rastr Whether the plot should be rasterized
#' @param quantCut If this is not NULL, a quantile cut is performed to threshold
#' the top and bottom of the distribution of numerical values.
#' This prevents skewed color scales caused by strong outliers.
#' The format of this should be c(x,y) where x is the lower threshold and y is
#' the upper threshold. For example, quantileCut = c(0.025,0.975) will take the
#' 2.5th percentile and 97.5 percentile of values and set values below/above to
#' the value of the 2.5th and 97.5th percentile values respectively.
#' @param continuousSet The name of a continuous palette from ArchRPalettes
#' for visualizing colorBy in the embedding if a continuous color set is desired.
#' @param discreteSet The name of a discrete palette from ArchRPalettes for
#' visualizing colorBy in the embedding if a discrete color set is desired.
#' @param randomize A boolean value that indicates whether to randomize points
#' prior to plotting to prevent cells from one cluster being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y
#' axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param addArrow A boolean value that indicates whether to add a smoothed
#' arrow in the embedding based on the aligned trajectory.
#' @param smoothWindow An integer value indicating the
#' smoothing window for creating inferred Arrow overlay on to embedding.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- TrajectoryPlot(
#'    object = obj,
#'    reduction = "dm",
#'    continuousSet = "blueYellow",
#'    size = 1
#' )
#' }
TrajectoryPlot <- function(object = NULL,
                           trajectory = "Trajectory",
                           reduction = NULL,
                           #name = "Trajectory",
                           size = 0.2,
                           rastr = FALSE,
                           quantCut = c(0.01, 0.99),
                           #quantHex = 0.5,
                           continuousSet = NULL,
                           discreteSet = NULL,
                           randomize = TRUE,
                           keepAxis = FALSE,
                           baseSize = 6,
                           addArrow = FALSE,
                           smoothWindow = 5) {
  dfT <- DataFrame(object@meta.data[, trajectory])
  rownames(dfT) <- colnames(object)
  idxRemove <- which(is.na(dfT[, 1]))
  df <-
    DataFrame(object@reductions[[reduction]]@cell.embeddings[, 1:2])
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("DM1", "DM2", "PseudoTime")


  plotParams <- list()
  plotParams$x <- df[, 1]
  plotParams$y <- df[, 2]
  plotParams$title <- reduction
  plotParams$baseSize <- baseSize
  plotParams$color <- as.vector(dfT$PseudoTime)
  plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
  plotParams$continuousSet <- "horizonExtra"
  plotParams$discreteSet <- "stallion"
  plotParams$title <- trajectory
  plotAs <- "points"
#    paste(plotParams$title, " colored by\ncolData : ", name)
#  if (is.null(plotAs)) {
#    plotAs <- "hexplot"
#  }

  plotParams$xlabel <- "DM1"
  plotParams$ylabel <- "DM2"

  if (!is.null(continuousSet)) {
    plotParams$continuousSet <- continuousSet
  }
  if (!is.null(continuousSet)) {
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  plotParams$color <- as.vector(plotParams$color)


  out <- do.call(ArchR::ggPoint, plotParams)
  out <- out +
    theme_cowplot() + ggtitle(trajectory) +
    xlab(colnames(df)[1]) + ylab(colnames(df)[2]) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "cm"),
    )


  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  colnames(dfT) <- c("x", "y", "PseudoTime", "value")
  dfT <- as.data.frame(dfT)

  if (addArrow) {
    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
      lapply(colMeans) %>%
      Reduce("rbind", .) %>%
      data.frame

    dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), , drop = FALSE])

    out <- out + geom_path(
      data = dfArrow,
      aes(x, y, color = NULL),
      size = 1,
      arrow = arrow(type = "open", length = unit(0.1, "inches"))
    )
  }

  return(out)
}


#' Plot trajectory heatmap
#'
#' This function generates a heatmap to visualize the data along inferred trajectory.
#' It was modified from the package ArchR \code{\link{plotTrajectoryHeatmap}}.
#' For more details, check here \url{https://www.archrproject.com/reference/plotTrajectoryHeatmap.html}.
#'
#'
#' @param trajectory A SummarizedExperiment object that results from calling \code{\link{GetTrajectory}}
#' @param varCutOff The "Variance Quantile Cutoff" to be used for identifying
#' the top variable features across the given trajectory.
#' Only features with a variance above the provided quantile will be retained.
#' @param maxFeatures The maximum number of features, ordered by variance,
#' to consider from useMatrix when generating a trajectory
#' @param scaleRows A boolean value that indicates whether row-wise z-scores
#' should be computed on the matrix
#' @param rowOrder If wanting to set the order of rows to be plotted, the
#' indices (integer or character correpsonding to rownmaes) can be provided here.
#' @param limits A numeric vector of two numbers that represent the lower and
#' upper limits of the heatmap color scheme.
#' @param labelRows A boolean value that indicates whether all rows should be
#' labeled on the side of the heatmap.
#' @param pal A custom continuous palette used to override the default
#' continuous palette for the heatmap.
#' @param labelMarkers A character vector listing the rownames that should be
#' labeled on the side of the heatmap.
#' @param labelTop A number indicating how many of the top N features,
#' based on variance of the matrix
#' @param name Name of the matrix
#' @param returnMatrix A boolean value that indicates whether the final heatmap
#' matrix should be returned in lieu of plotting the actual heatmap.
#'
#' @return a heatmap
#' @export
#'
TrajectoryHeatmap <- function(trajectory,
                              varCutOff = 0.9,
                              maxFeatures = 25000,
                              scaleRows = TRUE,
                              rowOrder = NULL,
                              limits = c(-1.5, 1.5),
                              labelRows = FALSE,
                              pal = NULL,
                              labelMarkers = NULL,
                              labelTop = 50,
                              name = "Heatmap",
                              returnMatrix = FALSE) {
  mat <- assay(trajectory)

  #Rows with NA
  rSNA <- rowSums(is.na(mat))
  if (sum(rSNA > 0) > 0) {
    message("Removing rows with NA values...")
    mat <- mat[rSNA == 0, ]#Remove NA Rows
  }

  varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))

  orderedVar <- FALSE

  if (is.null(rowOrder)) {
    mat <- mat[order(varQ, decreasing = TRUE), ]
    orderedVar <- TRUE
    if (is.null(varCutOff) & is.null(maxFeatures)) {
      n <- nrow(mat)
    } else if (is.null(varCutOff)) {
      n <- maxFeatures
    } else if (is.null(maxFeatures)) {
      n <- (1 - varCutOff) * nrow(mat)
    } else{
      n <- min((1 - varCutOff) * nrow(mat), maxFeatures)
    }
    n <- min(n, nrow(mat))
    mat <- mat[head(seq_len(nrow(mat)), n),]
  }


  if (!is.null(labelTop)) {
    if (orderedVar) {
      idxLabel <- rownames(mat)[seq_len(labelTop)]
    } else{
      idxLabel <-
        rownames(mat)[order(varQ, decreasing = TRUE)][seq_len(labelTop)]
    }
  } else{
    idxLabel <- NULL
  }

  if (scaleRows) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }
  if (nrow(mat) == 0) {
    stop("No Features Remaining!")
  }

  if (is.null(pal)) {
    pal <- ArchR::paletteContinuous(set = "blueYellow", n = 100)
  }

  if (!is.null(rowOrder)) {
    idx <- rowOrder
  } else{
    idx <- order(apply(mat, 1, which.max))
  }

  ht <- ArchR:::.ArchRHeatmap(
    mat = mat[idx, ],
    scale = FALSE,
    limits = c(min(mat), max(mat)),
    color = pal,
    clusterCols = FALSE,
    clusterRows = FALSE,
    labelRows = labelRows,
    labelCols = FALSE,
    customRowLabel = match(idxLabel, rownames(mat[idx,])),
    showColDendrogram = TRUE,
    name = name,
    draw = FALSE
  )

  if (returnMatrix) {
    return(mat[idx, ])
  } else{
    return(ht)
  }
}
