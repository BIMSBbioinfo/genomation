\docType{methods}
\name{scoreMatrixBin}
\alias{scoreMatrixBin}
\alias{scoreMatrixBin,GRanges,GRanges-method}
\alias{scoreMatrixBin,RleList,GRanges-method}
\alias{scoreMatrixBin,modRleList,GRanges-method}
\title{get scores overlapping with windows in a scoreMatrix object}
\usage{
  scoreMatrixBin(target,windows,bin.num=10,bin.op="mean",strand.aware=FALSE,...)
}
\arguments{
  \item{target}{a \code{RleList} or a \code{modRleList} or
  \code{GRanges} object to be overlapped with ranges in
  \code{windows}}

  \item{windows}{a \code{GRanges} object that contains the
  windows of interest. It could be promoters, CpG islands,
  exons, introns. However the sizes of windows does NOT
  have to be equal.}

  \item{bin.num}{A single \code{integer} value denoting how
  many bins there should be for each window}

  \item{bin.op}{A bin operation that is either one of the
  following strings: "max","min","mean". The operation is
  applied on the values in the bin. Defaults to "mean"}

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
  \code{windows} can be a predefined region such as CpG
  islands or gene bodies that are not necessarily
  equi-width. Each window will be chopped to equal number
  of bins based on \code{bin.num} option.
}
\seealso{
  \code{\link{scoreMatrix}},\code{\link{modCoverage}}
}

