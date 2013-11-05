\docType{methods}
\name{ScoreMatrixList}
\alias{ScoreMatrixList}
\title{ScoreMatrixList constructor}
\usage{
  ScoreMatrixList(l, granges=NULL, bin.num=NULL,
    bin.op='mean', strand.aware=FALSE, ...)
}
\arguments{
  \item{l}{can be a list of \code{scoreMatrix} objects,
  that are coerced to the \code{ScoreMatrixList}, a list of
  \code{RleList} objects, or a character vector specifying
  the locations of mulitple bam files that are used to
  construct the \code{scoreMatrixList}. If l is either a
  RleList object or a character vector of files, it is
  obligatory to give a granges argument.}

  \item{granges}{a \code{GenomicRanges} containing
  viewpoints for the scoreMatrix or ScoreMatrixList
  functions}

  \item{bin.num}{an integer telling the number of bins to
  bin the score matrix}

  \item{bin.op}{an name of the function that will be used
  for smoothing windows of ranges}

  \item{strand.aware}{a boolean telling the function
  whether to reverse the coverage of ranges that come from
  - strand (e.g. when plotting enrichment around
  transcription start sites)}

  \item{...}{other arguments that can be passed to the
  functio}
}
\value{
  returns a \code{ScoreMatrixList} object
}
\description{
  Construct a list of scoreMatrixObjects that can be used
  for plotting
}

