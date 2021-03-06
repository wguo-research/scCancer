% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAnnotation.R
\name{runCellInteraction}
\alias{runCellInteraction}
\title{runCellInteraction}
\usage{
runCellInteraction(
  expr,
  cellSetName = "default",
  species = "human",
  savePath = NULL
)
}
\arguments{
\item{expr}{A Seurat object.}

\item{cellSetName}{The colunm name of `expr`'s `meta.data`, used to indicate the cell set annotation.}

\item{species}{A character string indicating what species the sample belong to.
Only "human"(default) or "mouse" are allowed.}

\item{savePath}{A path to save the results files. If NULL, the 'statPath' will be used instead.}
}
\value{
A data frame which contains the cell sets ligand-receptor pairs and their scores.
}
\description{
runCellInteraction
}
