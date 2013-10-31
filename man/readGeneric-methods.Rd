\docType{methods}
\name{readGeneric}
\alias{readGeneric}
\alias{readGeneric,character-method}
\title{Read a tabular file and convert it to GRanges.}
\usage{
  readGeneric(location,remove.unsual=T)
}
\arguments{
  \item{file}{location of the file, a character string such
  as: "/home/user/my.bed"}

  \item{header}{whether the original file contains a header
  line which designates the column names.}

  \item{col.names}{a named \code{list} that maps individual
  columns to chr, start, end, strand information. e.g.
  list(chr=4, start=6, end=7, strand=9, score=10). If
  header = TRUE, colnames parameter does not have any
  effect. If the keep.medata is set to FALSE, only
  chromosome, start, end and strand colums will be used to
  construct the GRanges object, otherwise, all columns will
  be read in.}

  \item{starts.in.df.are.0based}{a boolean which tells
  whether the ranges in the bed file are 0 or 1 base
  encoded. A 0 based encoding is persumed}

  \item{remove.unsual}{if TRUE(default) remove the
  chromomesomes with unsual names, such as chrX_random}

  \item{sep}{a single character which designates the
  separator in the file. The default value is tab.}
}
\value{
  \code{\link{GRanges}} object
}
\description{
  The function can take a bed file with or without a
  predesignated header, or any other tabular file. If the
  file is in bed format, the user must specifiy which
  columns contain chromosome, start, and end information,
  using the col.names argument. Strand information is not
  compulsory.
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}

