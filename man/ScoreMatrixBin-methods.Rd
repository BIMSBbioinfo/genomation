\docType{methods}
\name{ScoreMatrixBin}
\alias{ScoreMatrixBin}
\alias{ScoreMatrixBin,character,GRanges-method}
\alias{ScoreMatrixBin,GRanges,GRanges-method}
\alias{ScoreMatrixBin,RleList,GRanges-method}
\title{Get bin score for bins on each window}
\usage{
  ScoreMatrixBin(target, windows, bin.num = 10,
    bin.op = "mean", strand.aware = FALSE,
    weight.col = NULL, is.noCovNA = FALSE, ...)
}
\arguments{
  \item{target}{\code{RleList} or a \code{modRleList} or
  \code{GRanges} object to be overlapped with ranges in
  \code{windows}}

  \item{windows}{\code{GRanges} object that contains the
  windows of interest.  It could be promoters, CpG islands,
  exons, introns. However the sizes of windows does NOT
  have to be equal.}

  \item{bin.num}{single \code{integer} value denoting how
  many bins there should be for each window}

  \item{bin.op}{bin operation that is either one of the
  following strings: "max","min","mean". The operation is
  applied on the values in the bin. Defaults to "mean"}

  \item{strand.aware}{If TRUE (default: FALSE), the strands
  of the windows will be taken into account in the
  resulting \code{scoreMatrix}.  If the strand of a window
  is -, the values of the bins for that window will be
  reversed}

  \item{weight.col}{if the object is \code{GRanges} object
  a numeric column in meta data part can be used as
  weights. This is particularly useful when genomic regions
  have scores other than their coverage values, such as
  percent methylation, conservation scores, GC content,
  etc.}

  \item{is.noCovNA}{(Default:FALSE) if TRUE,and if 'target'
  is a GRanges object with 'weight.col' provided, the bases
  that are uncovered will be preserved as NA in the
  returned object. This useful for situations where you can
  not have coverage all over the genome, such as CpG
  methylation values.}
}
\value{
  returns a \code{scoreMatrix} object
}
\description{
  The function firsts bins each window to equal number of
  bins, and calculates the a summary metrix for scores of
  each bin (currently, mean, max and min supported) A
  scoreMatrix object can be used to draw average profiles
  or heatmap of read coverage or wig track-like data.
  \code{windows} can be a predefined region such as CpG
  islands or gene bodies that are not necessarily
  equi-width. Each window will be chopped to equal number
  of bins based on \code{bin.num} option.
}
\examples{
data(cage)
  data(cpgi)
  myMat=ScoreMatrixBin(target=cage,
                       windows=cpgi,bin.num=10,bin.op="mean",weight.col="tpm")
  plot(colMeans(myMat,na.rm=TRUE),type="l")

  myMat2=ScoreMatrixBin(target=cage,
                        windows=promoters,bin.num=10,bin.op="mean",
                        weight.col="tpm",strand.aware=TRUE)
  plot(colMeans(myMat2,na.rm=TRUE),type="l")
}
\seealso{
  \code{\link{ScoreMatrix}}
}

