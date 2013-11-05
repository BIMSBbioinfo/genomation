\docType{methods}
\name{scaleScoreMatrixList}
\alias{scaleScoreMatrixList}
\title{Scale the ScoreMatrixList}
\usage{
  scaleScoreMatrixList(sml, columns, rows, scalefun, ...)
}
\arguments{
  \item{sml}{a \code{ScoreMatrixList} object}

  \item{columns}{a \code{columns} whether to scale the
  matrix by columns. Set by default to FALSE}

  \item{rows}{a \code{rows} Whether to scale the matrix by
  rows. Set by default to TRUE}

  \item{scalefun}{a function object that takes as input a
  matrix and returns a matrix.}

  \item{...}{other argments that be passed to the function
  By default the argument is set to the R scale function
  with center=TRUE and scale=TRUE}
}
\value{
  \code{ScoreMatrixList} object
}
\description{
  Scales each scoreMatrix in the ScoreMatrixList object, by
  rows and/or columns
}

