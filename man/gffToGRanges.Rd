\docType{methods}
\name{gffToGRanges}
\alias{gffToGRanges}
\title{Converts a gff formated data.frame into a GenomicRanges object.
The GenomicRanges object needs to be properly formated for the function to work.}
\usage{
  gffToGRanges(gff.file, split.group = FALSE,
    split.char = ";", filter = NULL, zero.based = FALSE)
}
\arguments{
  \item{gff.file}{path to a gff formatted file}

  \item{split.group}{boolean, whether to split the 9th
  column of the file}

  \item{split.char}{character that is used as a separator
  of the 9th column. ';' by default}

  \item{filter}{a character designating which elements to
  retain from the gff file (e.g. exon, CDS, ...)}

  \item{zero.based}{\code{boolean} whether the coordinates
  are 0 or 1 based. 0 is the default}
}
\value{
  returns a \code{GenomicRanges} object
}
\description{
  Converts a gff formated data.frame into a GenomicRanges
  object. The GenomicRanges object needs to be properly
  formated for the function to work.
}
\examples{
gff.file = system.file('extdata/chr21.refseq.hg19.gtf', package='genomation')
gff = gffToGRanges(gff.file, split.group=TRUE)
}

