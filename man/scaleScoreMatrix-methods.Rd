\docType{methods}
\name{scaleScoreMatrix}
\alias{scaleScoreMatrix}
\alias{scaleScoreMatrix,ScoreMatrix-method}
\title{Scales the values in the matrix by rows and/or columns}
\usage{
  scaleScoreMatrix(mat, columns=FALSE, rows=TRUE,
    scalefun=NULL, ...)

  \S4method{scaleScoreMatrix}{ScoreMatrix}(mat,
    columns = FALSE, rows = TRUE, scalefun = NULL, ...)
}
\arguments{
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
}
\description{
  Scales the values in the matrix by rows and/or columns
}

