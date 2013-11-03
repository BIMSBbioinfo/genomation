\docType{methods}
\name{ScoreMatrixList}
\alias{order}
\alias{scaleScoreMatrixList}
\alias{ScoreMatrixList}
\alias{unionScoreMatrixList}
\title{ScoreMatrixList constructor}
\usage{
  ScoreMatrixList(l, granges, bin)
}
\arguments{
  \item{l:}{can be a list of \code{scoreMatrix} objects,
  that are coerced to the \code{ScoreMatrixList}, a list of
  \code{RleList} objects, or a character vector specifying
  the locations of mulitple bam files that are used to
  construct the \code{scoreMatrixList}. If l is either a
  RleList object or a character vector of files, it is
  obligatory to give a granges argument.}

  \item{granges:}{a \code{GenomicRanges} containing
  viewpoints for the scoreMatrix or ScoreMatrixList
  functions}

  \item{bin:}{an integer telling the number of bins to bin
  the score matrix}

  \item{sml}{a \code{ScoreMatrixList} object}

  \item{columns}{a \code{columns} whether to scale the
  matrix by columns. Set by default to FALSE.}

  \item{rows}{a \code{rows} Whether to scale the matrix by
  rows. Set by default to TRUE}

  \item{scalefun}{a function object that takes as input a
  matrix and returns a matrix.  By default the argument is
  set to the R scale function with center=TRUE and
  scale=TRUE}

  \item{sml}{a \code{ScoreMatrixList} object}

  \item{sml}{a \code{ScoreMatrixList} object}

  \item{ord.vec}{an integer vector}
}
\value{
  returns a \code{ScoreMatrixList} object

  \code{ScoreMatrixList} object

  \code{ScoreMatrixList} object

  \code(ScoreMatrixList) object
}
\description{
  Construct a list of scoreMatrixObjects that can be used
  for plotting

  Scales each scoreMatrix in the ScoreMatrixList object

  Returns a union of rows for each matrix in a
  ScoreMatrixList object. This is done using the rownames
  of each element in the list.

  Reorder all elements of a ScoreMatrixList to a given
  ordering vector
}

