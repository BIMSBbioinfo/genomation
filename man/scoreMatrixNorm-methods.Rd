\docType{methods}
\name{scoreMatrixNorm}
\alias{scoreMatrixNorm}
\alias{scoreMatrixNorm,GRanges,GRanges,ANY-method}
\alias{scoreMatrixNorm,GRanges,RleList-method}
\title{get scores overlapping with windows in a scoreMatrix object}
\usage{
  scoreMatrix(target,windows)
}
\arguments{
  \item{target}{a \code{RleList} or a \code{modRleList}
  object to be overlapped with ranges in \code{windows}}

  \item{windows}{a \code{GRanges} object that will be
  randomly placed across the genome and overlap of these
  random regions with \code{target} will be the background
  distribution of association between \code{target} and
  \code{query}.}

  \item{bin.num}{A single \code{integer} value denoting how
  many bins there should be for each window}

  \item{...}{parameters to be passed to \code{modCoverage}
  function}
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

