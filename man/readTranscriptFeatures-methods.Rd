\docType{methods}
\name{readTranscriptFeatures}
\alias{readTranscriptFeatures}
\alias{readTranscriptFeatures,character-method}
\title{Function for reading exon intron and promoter structure from a given bed file}
\usage{
readTranscriptFeatures(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE)

\S4method{readTranscriptFeatures}{character}(location, remove.unsual, up.flank,
  down.flank, unique.prom)
}
\arguments{
  \item{location}{location of the bed file with 12 or more
  columns}

  \item{remove.unsual}{remove the chromomesomes with unsual
  names, mainly random chromsomes etc}

  \item{up.flank}{up-stream from TSS to detect promoter
  boundaries}

  \item{down.flank}{down-stream from TSS to detect promoter
  boundaries}

  \item{unique.prom}{get only the unique promoters,
  promoter boundaries will not have a gene name if you set
  this option to be TRUE}
}
\value{
a \code{\link{GRangesList}} containing locations of
exon/intron/promoter/TSS
}
\description{
Function for reading exon intron and promoter structure
from a given bed file
}
\note{
one bed track per file is only accepted, the bed files with
multiple tracks will cause en error
}
\examples{
my.bed12.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
  my.bed12.file
  feats = readTranscriptFeatures(my.bed12.file)
  names(feats)
  sapply(feats, head)
}

