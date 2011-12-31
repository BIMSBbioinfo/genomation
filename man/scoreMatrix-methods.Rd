\docType{methods}
\name{scoreMatrix}
\alias{scoreMatrix}
\alias{scoreMatrix,GRanges,GRanges,ANY-method}
\alias{scoreMatrix,GRanges,RleList-method}
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

  \item{...}{parameters to be passed to \code{modCoverage}
  function}
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

