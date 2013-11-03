\docType{methods}
\name{gffToGRanges}
\alias{gffToGRanges}
\title{Converts a gff formated data.frame into a GenomicRanges object.
The GenomicRanges object needs to be properly formated for the function to work. If some of the ranges have negative width, then the function will issue a warning (happens sometimes with ranges downloaded from Ensembl)}
\usage{
  gffToGRanges(gff)
}
\arguments{
  \item{gff:}{a gff3 formated \code{data.frame}}

  \item{filter:}{a character designating which elements to
  retain from the gff file (e.g. exon, CDS, ...)}
}
\value{
  returns a \code{GenomicRanges} object
}
\description{
  Converts a gff formated data.frame into a GenomicRanges
  object. The GenomicRanges object needs to be properly
  formated for the function to work. If some of the ranges
  have negative width, then the function will issue a
  warning (happens sometimes with ranges downloaded from
  Ensembl)
}

