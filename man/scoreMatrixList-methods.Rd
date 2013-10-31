\docType{methods}
\name{ScoreMatrixList}
\alias{ScoreMatrixList}
\title{ScoreMatrixList constructor}
\usage{
  ScoreMatrixList(l, granges = NULL, bin = NULL, ...)
}
\arguments{
  \item{l}{corresponds to can be a list of
  \code{scoreMatrix} objects, that are coerced to the
  \code{ScoreMatrixList}, or a list of \code{RleList}
  objects that are used to construct the
  \code{scoreMatrixList}}

  \item{granges}{a \code{GenomicRanges} containing
  viewpoints for the scoreMatrix or ScoreMatrixList
  functions}

  \item{bin}{an integer telling the number of bins to bin
  the score matrix}
}
\value{
  returns a \code{ScoreMatrixList} object
}
\description{
  Construct a list of scoreMatrixObjects that can be used
  for plotting
}

