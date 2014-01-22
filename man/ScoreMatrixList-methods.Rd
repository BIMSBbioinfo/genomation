\docType{methods}
\name{ScoreMatrixList}
\alias{ScoreMatrixList}
\title{Make ScoreMatrixList from multiple targets}
\usage{
ScoreMatrixList(targets, windows = NULL, bin.num = NULL, bin.op = "mean",
  strand.aware = FALSE, weight.col = NULL, is.noCovNA = FALSE,
  type = "", rpm = FALSE, unique = FALSE, extend = 0, param = NULL)
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

  \item{weight.col}{if the object is \code{GRanges} object
  a numeric column in meta data part can be used as
  weights. This is particularly useful when genomic regions
  have scores other than their coverage values, such as
  percent methylation, conservation scores, GC content,
  etc.}

  \item{is.noCovNA}{(Default:FALSE) if TRUE,and if
  'targets' is a GRanges object with 'weight.col' provided,
  the bases that are uncovered will be preserved as NA in
  the returned object. This useful for situations where you
  can not have coverage all over the genome, such as CpG
  methylation values.}

  \item{type}{if \code{targets} is a character vector of
  file paths, then type designates the type of the
  corresponding files (bam or bigWig)}

  \item{rpm}{boolean telling whether to normalize the
  coverage to per milion reads.  FALSE by default.}

  \item{unique}{boolean which tells the function to remove
  duplicated reads based on chr, start, end and strand}

  \item{extend}{numeric which tells the function to extend
  the features ( i.e aligned reads) to total length
  ofwidth+extend}

  \item{param}{ScanBamParam object}
}
\value{
returns a \code{ScoreMatrixList} object
}
\description{
The function constructs a list of \code{ScoreMatrix}
objects in the form of \code{ScoreMatrixList} object. This
object can be visualized using \code{multiHeatMatrix},
\code{heatMeta} or \code{plotMeta}
}
\examples{
# visualize the distribution of cage clusters and cpg islands around promoters
 library(GenomicRanges)
 data(cage)
 data(cpgi)
 data(promoters)

 cage$tpm = NULL
 targets = GRangesList(cage=cage, cpgi=cpgi)
 sml = ScoreMatrixList(targets, promoters, bin.num=10, strand.aware=TRUE)
 multiHeatMatrix(sml)
}

