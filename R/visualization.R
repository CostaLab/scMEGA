#' Plot dynamic changes along pseudotime trajectory
#'
#' This function will generate a line plot to visualize the TF binding activity,
#' TF expression, and TF target expression along the trajectory.
#'
#' @param object A Seurat object used as input.
#' @param tf.use Which TF to plot.
#' @param tf.assay Assay name for TF binding activity. Default: "chromvar"
#' @param rna.assay Assay name for TF expression. Default: "RNA"
#' @param atac.assay Assay name for Peaks. Default: "ATAC"
#' @param target.assay Assay name for TF target expression. Default: "target"
#' @param trajectory.name Trajectory name for visualization.
#' @param groupEvery The number of sequential percentiles to group together when generating a trajectory.
#' This is similar to smoothing via a non-overlapping sliding window across pseudo-time.
#'
#' @return A ggplot object
#' @export
#'
PseudotimePlot <- function(object, tf.use,
                           tf.assay="chromvar",
                           rna.assay = "RNA",
                           atac.assay = "ATAC",
                           target.assay = "target",
                           trajectory.name = "Trajectory",
                           groupEvery=1
                           ){

    trajMM <- GetTrajectory(
        object,
        assay = tf.assay,
        trajectory.name=trajectory.name,
        groupEvery = groupEvery,
        slot = "data",
        smoothWindow = 7,
        log2Norm = FALSE
    )

    rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names

    df.tf.activity <- assay(trajMM)
    df.tf.activity <- t(scale(t(df.tf.activity)))
    df.tf.activity <- as.data.frame(df.tf.activity)
    colnames(df.tf.activity) <- seq(0, 100, groupEvery)[-1]
    df.tf.activity$tf <- toupper(rownames(df.tf.activity))
    df.tf.activity <- tidyr::pivot_longer(df.tf.activity, -tf,
                                              names_to = "pseudotime",
                                              values_to = "value")
    df.tf.activity$pseudotime <- as.numeric(df.tf.activity$pseudotime)

    trajGEX <- GetTrajectory(
        object,
        assay = rna.assay,
        trajectory.name =trajectory.name,
        groupEvery = groupEvery,
        slot = "data",
        smoothWindow = 7,
        log2Norm = TRUE
    )

    df.tf.expression <- assay(trajGEX)
    df.tf.expression <- t(scale(t(df.tf.expression)))
    df.tf.expression <- as.data.frame(df.tf.expression)
    colnames(df.tf.expression) <- seq(0, 100, groupEvery)[-1]
    df.tf.expression$tf <- toupper(rownames(df.tf.expression))
    df.tf.expression <- tidyr::pivot_longer(df.tf.expression, -tf,
                                              names_to = "pseudotime",
                                              values_to = "value")
    df.tf.expression$pseudotime <- as.numeric(df.tf.expression$pseudotime)

    traj.target <- GetTrajectory(
        object,
        assay = target.assay,
        trajectory.name=trajectory.name,
        groupEvery = groupEvery,
        slot = "data",
        smoothWindow = 7,
        log2Norm = FALSE
    )

    df.target <- assay(traj.target)
    df.target <- t(scale(t(df.target)))
    df.target <- as.data.frame(df.target)
    colnames(df.target) <- seq(0, 100, groupEvery)[-1]
    df.target$tf <- toupper(rownames(df.target))
    df.target <- tidyr::pivot_longer(df.target, -tf,
                                              names_to = "pseudotime",
                                              values_to = "value")
    df.target$pseudotime <- as.numeric(df.target$pseudotime)

    df.tf.activity$data <- "TF activity"
    df.tf.expression$data <- "TF expression"
    df.target$data <- "Targets expression"

    df.tf <- rbind(df.tf.activity, df.tf.expression, df.target)
    df.plot <- subset(df.tf, tf == tf.use)

    p <- ggplot(df.plot, aes(x = pseudotime, y = value, color = data)) +
         geom_smooth(method = "loess",se = FALSE) +
        ggtitle(tf.use) +
        cowplot::theme_cowplot() +
        ylab("") +
        theme(legend.title = element_blank())

    return(p)

}

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

  counts <- NULL
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
#' @param method.test Name of method for statistical test. Default: "wilcox.test"
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats as.formula
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

    counts <- NULL
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
        geom_boxplot(aes_string(color = sample.annotation), outlier.shape = NA) +
        geom_jitter(aes_string(color = sample.annotation), size=0.4, alpha=1) +
        facet_wrap(as.formula(paste("~", group.by)), nrow = 1) +
        ggpubr::stat_compare_means(comparisons = comparisons,
                                   method = method.test) +
        cowplot::theme_cowplot() +
        xlab("") + ylab("") +
        theme(axis.text.x = element_blank(),
              legend.title = element_blank())

    } else{
      p <-
        ggplot(data = df, aes_string(x = sample.annotation, y = proportion)) +
        geom_boxplot(aes_string(color = sample.annotation), outlier.shape = NA) +
        geom_jitter(aes_string(color = sample.annotation), size=0.4, alpha=1) +
        facet_wrap(as.formula(paste("~", group.by)), nrow = 1) +
        cowplot::theme_cowplot() +
        xlab("") + ylab("") +
        theme(axis.text.x = element_blank(),
              legend.title = element_blank())
    }

    return(p)

  }


#' Plot correlation
#'
#' This function will generate a scatter plot to visualize the correlation between
#' expression and activity for all TFs passing the threshold.
#'
#' @param df A data frame generated by the function \code{\link{GetCorrelation}}
#'
#' @return A ggplot object
#' @export
#'
CorrelationPlot <- function(df) {
  p <-
    ggplot(data = df, aes(x = reorder(tfs,-correlation), y = correlation)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = tfs),
                             max.overlaps = 100) +
    xlab("TFs") + ylab("Correlation") +
    ggtitle(glue::glue("Number of TFs: {nrow(df_cor)}")) +
    theme_cowplot() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

  return(p)

}

#' Plot correlation heatmap
#'
#' This function will generate a combined heatmap with each one corresponding to
#' the TF expression and binding activity as selected by using the function
#' \code{\link{GetCorrelation}}. It was modified from the package \code{ArchR}.
#' More details can be found:
#' \url{https://www.archrproject.com/bookdown/myeloid-trajectory-monocyte-differentiation.html}
#'
#' @param trajectory1 A SummarizedExperiment object for TF activity
#' generated by the function \code{\link{GetTrajectory}}
#' @param trajectory2 A SummarizedExperiment object for TF expression
#' @param name1 Name for the first object
#' @param name2 Name for the second object
#' @param labelRows1 Whether or not to show all labels for the first heatmap
#' @param labelRows2 Whether or not to show all labels for the second heatmap
#' @param labelTop1 A number indicating how many of the top N features,
#' based on variance, should be labeled on the side of the first heatmap
#' @param labelTop2 A number indicating how many of the top N features,
#' based on variance, should be labeled on the side of the second heatmap
#' @param limits1 A numeric vector of two numbers that represent the lower and
#' upper limits of the second heatmap color scheme.
#' generated by the function \code{\link{GetTrajectory}}
#' @param limits2 A numeric vector of two numbers that represent the lower and
#' upper limits of the second heatmap color scheme.
#' generated by the function \code{\link{GetTrajectory}}
#' @return A heatmap
#' @export
#'
CorrelationHeatmap <- function(trajectory1,
                               trajectory2,
                               name1 = NULL,
                               name2 = NULL,
                               labelRows1 = TRUE,
                               labelRows2 = TRUE,
                               labelTop1 = 50,
                               labelTop2 = 50,
                               limits1 = c(-2, 2),
                               limits2 = c(-2, 2)) {
  trajCombined <- trajectory1

  assay(trajCombined, withDimnames = FALSE) <-
    t(apply(assay(trajectory2), 1, scale)) +
    t(apply(assay(trajectory1), 1, scale))

  combinedMat <- TrajectoryHeatmap(trajCombined,
                                   returnMatrix = TRUE,
                                   varCutOff = 0)

  rowOrder <- match(rownames(combinedMat), rownames(trajectory1))

  ht1 <- TrajectoryHeatmap(
    trajectory1,
    pal = paletteContinuous(set = "solarExtra"),
    varCutOff = 0,
    maxFeatures = nrow(trajectory1),
    rowOrder = rowOrder,
    limits = limits1,
    labelRows = labelRows1,
    labelTop = labelTop1,
    name = name1
  )

  ht2 <- TrajectoryHeatmap(
    trajectory2,
    pal = paletteContinuous(set = "horizonExtra"),
    varCutOff = 0,
    maxFeatures = nrow(trajectory2),
    rowOrder = rowOrder,
    limits = limits2,
    labelRows = labelRows2,
    labelTop = labelTop2,
    name = name2
  )

  ht <- ht1 + ht2

  return(ht)

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
#' @import dplyr
#' @import ggplot2
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
  dfT <- object@meta.data[, trajectory] %>%
    as.data.frame()

  rownames(dfT) <- colnames(object)
  idxRemove <- which(is.na(dfT[, 1]))
  df <- object@reductions[[reduction]]@cell.embeddings[, 1:2] %>%
    as.data.frame()

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
    cowplot::theme_cowplot() +
    ggtitle(trajectory) +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2]) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "cm")
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
#' @importFrom utils head
#' @importFrom SummarizedExperiment assay
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
  #colum all zero
  rSZero <- which(matrixStats::rowSds(mat) != 0)
  if(length(rSZero) < nrow(mat)){
    message("Removing rows without peaks...")
    mat <- mat[rSZero, ]#Remove 0 Rows
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


  if (!is.null(labelTop) & labelTop > 0) {
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

  if(!is.null(idxLabel)){
    customRowLabel <- match(idxLabel, rownames(mat[idx,]))
  } else{
      customRowLabel <- NULL
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
    customRowLabel = customRowLabel,
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
