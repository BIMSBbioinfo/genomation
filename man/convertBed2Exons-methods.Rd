\docType{methods}
\name{convertBed2Exons}
\alias{convertBed2Exons}
\alias{convertBed2Exons,data.frame-method}
\title{convert a data frame read-in from a bed file to a GRanges object for exons}
\usage{
  convertBed2Exons(bed.df)

  \S4method{convertBed2Exons}{data.frame}(bed.df)
}
\arguments{
  \item{bed.df}{a data.frame where column order and content
  resembles a bed file with 12 columns}
}
\value{
  \code{\link{GRanges}} object
}
\description{
  convert a data frame read-in from a bed file to a GRanges
  object for exons
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}

