\docType{methods}
\name{binMatrix}
\alias{binMatrix}
\alias{binMatrix,ScoreMatrix-method}
\alias{binMatrix,ScoreMatrixList-method}
\title{Bins the columns of a matrix using a user provided function}
\usage{
  binMatrix(x, bin.num = NULL, fun = "mean", ...)

  \S4method{binMatrix}{ScoreMatrix}(x, bin.num = NULL,
    fun = "mean", ...)

  \S4method{binMatrix}{ScoreMatrixList}(x, bin.num = NULL,
    fun = "mean", ...)
}
\arguments{
  \item{x}{\code{ScoreMatrix} or a \code{ScoreMatrixList}
  object}

  \item{bin.num}{\code{integer} number of bins in the final
  matrix}

  \item{fun}{\code{character} vector or an anonymous
  function that will be used for binning}
}
\value{
  \code{ScoreMatrix} or \code{ScoreMatrixList} object
}
\description{
  Bins the columns of a matrix using a user provided
  function
}

