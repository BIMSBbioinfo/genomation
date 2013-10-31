\docType{methods}
\name{heatmapProfile}
\alias{heatmapProfile}
\alias{heatmapProfile,ScoreMatrixList-method}
\title{Plot functions for score matrix list}
\usage{
  heatmapProfile(mat.list, mat.cols=NULL, ...)
}
\arguments{
  \item{mat.list}{a \code{ScoreMatrixList} object}

  \item{mat.cols}{colors to be used for plotting}

  \item{xmarks}{an integer number to lable the thick marks
  on the x axis of each heatmap. By default it takes the
  values of -ncol/2, 0, ncol/2}

  \item{ymarks}{a vector of that will lable the thick marks
  on the y axis}

  \item{y.at}{a numeric vector that will specify the
  positions of the thick marks on the y axis}

  \item{xcex,}{ycex an integer number which controls the
  character expansion on x and y axis}

  \item{cex.main}{an integer number which controls the
  character expansion of the plot label}

  \item{mar}{a vector of length 5 which controls the size
  of the margins. The order is the following: below, left,
  up, right, spacing between consecutive plots}

  \item{use.names}{whether to use the names of the
  ScoreMatrixList object to label each plot}

  \item{main}{whether to use the names of the
  ScoreMatrixList object to label each plot}

  \item{xlab,}{ylab name to be used for the x/y axis}

  \item{...}{other options (obselete for now)}
}
\description{
  Functions plot a ScoreMatrixList object as a panel of
  heatmaps
}
\examples{
l = lapply(seq(20, 40,5), function(x)as(matrix(rpois(1000, x), ncol=25), 'scoreMatrix'))
}

