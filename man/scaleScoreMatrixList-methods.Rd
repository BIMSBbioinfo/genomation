\docType{methods}
\name{scaleScoreMatrixList}
\alias{scaleScoreMatrixList}
\title{Scale the ScoreMatrixList}
\usage{
  scaleScoreMatrixList(mat, columns=FALSE, rows=TRUE, ...)
}
\arguments{
  \item{sml}{a \code{ScoreMatrixList} object}

  \item{columns}{a \code{columns} whether to scale the
  matrix by columns. Set by default to FALSE.}

  \item{rows}{a \code{rows} Whether to scale the matrix by
  rows. Set by default to TRUE}

  \item{scalefun}{a function object that takes as input a
  matrix and returns a matrix.  By default the argument is
  set to the R scale function with center=TRUE and
  scale=TRUE}
}
\value{
  \code(ScoreMatrixList) object
}
\description{
  Scales each scoreMatrix in the ScoreMatrixList object
}

