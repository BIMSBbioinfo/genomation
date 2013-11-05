\docType{methods}
\name{binMatrix}
\alias{binMatrix}
\alias{binMatrix,ScoreMatrix-method}
\alias{scaleScoreMatrix}
\alias{scaleScoreMatrix,ScoreMatrix-method}
\title{Bins the columns of a matrix using a user provided function}
\usage{
  binMatrix(mat, nbins = NULL, fun = "mean", ...)

  \S4method{binMatrix}{ScoreMatrix}(mat, nbins = NULL,
    fun = "mean", ...)

  scaleScoreMatrix(mat, columns=FALSE, rows=TRUE,
    scalefun=NULL, ...)

  \S4method{scaleScoreMatrix}{ScoreMatrix}(mat,
    columns = FALSE, rows = TRUE, scalefun = NULL, ...)
}
\arguments{
  \item{mat}{\code{ScoreMatrix} object}

  \item{nbins}{\code{integer} number of bins in the final
  matrix}

  \item{fun}{\code{character} vector representing the
  function to be used for bining}

  \item{mat}{\code{ScoreMatrix} object}

  \item{columns}{\code{columns} whether to scale the matrix
  by columns. Set by default to FALSE.}

  \item{rows}{\code{rows} Whether to scale the matrix by
  rows. Set by default to TRUE}

  \item{scalefun}{function object that takes as input a
  matrix and returns a matrix. By default the argument is
  set to (x - mean(x))/(max(x)-min(x)+1)}
}
\value{
  \code{ScoreMatrix} object

  \code{ScoreMatrix} object
}
\description{
  Bins the columns of a matrix using a user provided
  function

  Scales the values in the matrix by rows and/or columns
}

