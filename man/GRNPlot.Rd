% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grn.R
\name{GRNPlot}
\alias{GRNPlot}
\title{Get a graph}
\usage{
GRNPlot(
  df.grn,
  tfs.use = NULL,
  show.tf.labels = TRUE,
  tfs.timepoint = NULL,
  genes.cluster = NULL,
  genes.use = NULL,
  genes.highlight = NULL,
  cols.highlight = "#984ea3",
  seed = 42,
  plot.importance = TRUE,
  min.importance = 2,
  remove.isolated = FALSE
)
}
\arguments{
\item{df.grn}{A data frame representing predicted network}

\item{tfs.timepoint}{Time points of TFs}

\item{genes.cluster}{A data frame containing clustering results of genes}

\item{genes.highlight}{A string list to include gene names for plotting}

\item{cols.highlight}{Color code for highlighted genes}

\item{seed}{Random seet}

\item{plot.importance}{Whether or not plot the scatter plot to visualize importance score of each TF}

\item{min.importance}{The minimum importance score for showing the TF labels.}
}
\value{
A ggplot object
}
\description{
This function will generate a graph to visualize the predicted gene regulatory network
}
