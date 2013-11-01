\docType{methods}
\name{ScoreMatrix}
\alias{scaleScoreMatrix}
\alias{ScoreMatrix}
\alias{ScoreMatrix,character,GRanges-method}
\alias{ScoreMatrix,GRanges,GRanges-method}
\alias{ScoreMatrix,RleList,GRanges-method}
\title{Get base-pair score for bases in each window}
\usage{
  ScoreMatrix(target, windows, strand.aware = FALSE,
    weight.col = NULL, is.noCovNA = FALSE)

  scaleScoreMatrix(mat, columns = FALSE, rows = TRUE,
    scalefun = function(x) scale(x), ...)
}
\arguments{
  \item{target}{a \code{RleList} , \code{GRanges} or a BAM
  file to be overlapped with ranges in \code{windows}}

  \item{windows}{a \code{GRanges} object that contains the
  windows of interest.  It could be promoters, CpG islands,
  exons, introns.  However the sizes of windows have to be
  equal.}

  \item{strand.aware}{If TRUE (default: FALSE), the strands
  of the windows will be taken into account in the
  resulting \code{ScoreMatrix}.  If the strand of a window
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

  \item{mat}{a \code{ScoreMatrix} object}

  \item{columns}{a \code{columns} whether to scale the
  matrix by columns. Set by default to FALSE.}

  \item{rows}{a \code{rows} Whether to scale the matrix by
  rows. Set by default to TRUE}

  \item{scalefun}{a function object that takes as input a
  matrix and returns a matrix. By default the argument is
  set to the R scale function with center=TRUE and
  scale=TRUE}
}
\value{
  returns a \code{ScoreMatrix} object

  \code{ScoreMatrix} object
}
\description{
  The funcion produces a base-pair resolution matrix of
  scores for given equal width windows of interest. The
  returned matrix can be used to draw meta profiles or
  heatmap of read coverage or wig track-like data. The
  \code{windows} argument can be a predefined region around
  transcription start sites or other regions of interest
  that have equal lengths The function removes all window
  that fall off the Rle object - have the start coordinate
  < 1 or end coordinate > length(Rle) The function takes
  the intersection of names in the Rle and GRanges objects

  Scales the values in the matrix by rows and/or columns
}
\examples{
# When target is GRanges
         data(cage)
         data(promoters)
         scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
                                 weight.col="tpm")

# When target is RleList
covs=coverage(cage)
scores2=ScoreMatrix(target=covs,windows=promoters,strand.aware=TRUE)


# When target is a bam file
 # bamfile="example.bam"
 # scores3=ScoreMatrix(target=bamfile,windows=promoters,strand.aware=TRUE)
}
\seealso{
  \code{\link{ScoreMatrixBin}}
}

