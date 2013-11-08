\docType{methods}
\name{ScoreMatrixList}
\alias{ScoreMatrixList}
\title{Make ScoreMatrixList from multiple targets}
\usage{
  ScoreMatrixList(targets, windows = NULL, bin.num = NULL,
    bin.op = "mean", strand.aware = FALSE, ...)
}
\arguments{
  \item{targets}{can be a list of \code{scoreMatrix}
  objects, that are coerced to the \code{ScoreMatrixList},
  a list of \code{RleList} objects, or a character vector
  specifying the locations of mulitple bam files that are
  used to construct the \code{scoreMatrixList}. If l is
  either a RleList object or a character vector of files,
  it is obligatory to give a granges argument.}

  \item{windows}{\code{GenomicRanges} containing viewpoints
  for the scoreMatrix or ScoreMatrixList functions}

  \item{bin.num}{an integer telling the number of bins to
  bin the score matrix}

  \item{bin.op}{an name of the function that will be used
  for smoothing windows of ranges}

  \item{strand.aware}{a boolean telling the function
  whether to reverse the coverage of ranges that come from
  - strand (e.g. when plotting enrichment around
  transcription start sites)}

  \item{type}{if \code{targets} is a character vector of
  file paths, then type designates the type of the
  corresponding files (bam or bigWig)}

  \item{...}{other arguments that can be passed to the
  function}
}
\value{
  returns a \code{ScoreMatrixList} object
}
\description{
  The function constructs a list of \code{ScoreMatrix}
  objects in the form of \code{ScoreMatrixList} object.
  This object can be visualized using
  \code{multiHeatMatrix}, \code{heatMeta} or
  \code{plotMeta}
}

