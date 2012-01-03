\docType{methods}
\name{scoreMatrix}
\alias{scoreMatrix}
\alias{scoreMatrix,GRanges,GRanges,ANY-method}
\alias{scoreMatrix,GRanges,RleList-method}
\title{get scores overlapping with windows in a scoreMatrix object}
\usage{
  scoreMatrix(target,windows,target,windows,strand.aware=FALSE,...)
}
\arguments{
  \item{target}{a \code{RleList} or a \code{modRleList} or
  \code{GRanges} object to be overlapped with ranges in
  \code{windows}}

  \item{windows}{a \code{GRanges} object that contains the
  windows of interest. It could be promoters, CpG islands,
  exons, introns. However the sizes of windows have to be
  equal.}

  \item{strand.aware}{If TRUE (default: FALSE), the strands
  of the windows will be taken into account in the
  resulting \code{scoreMatrix}. If the strand of a window
  is -, the values of the bins for that window will be
  reversed}

  \item{...}{parameters to be passed to \code{modCoverage}
  function. Only needed when target is \code{GRanges}.}
}
\value{
  returns a \code{scoreMatrix} object
}
\description{
  A scoreMatrix object can be used to draw average profiles
  or heatmap of read coverage or wig track-like data.
  \code{windows} can be a predefined region around
  transcription start sites or other regions of interest
  that have equal lengths
}
\seealso{
  \code{link{scoreMatrixBin-methods}},
  \code{link{modCoverage-methods}}
}

