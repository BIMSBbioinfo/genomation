\docType{methods}
\name{modCoverage}
\alias{modCoverage}
\alias{modCoverage,GRanges,character-method}
\title{gets a weighted coverage vector from a GRanges object with a numeric column}
\usage{
  modCoverage(x,col.name,multiply=1,add=0)
}
\arguments{
  \item{x}{GRanges object with elementData}

  \item{col.name}{name of the column that has a numeric
  value to be used as a value in RleList object. It should
  be a column from the \code{elementMetadata} section of
  the object.}

  \item{multiply}{a value (default:1000) that will be
  multiplied by the value denoted by \code{col.name} in
  GRanges object. This is useful to retain some of the
  values with decimals.}

  \item{add}{This value (default:1) be added to the value
  denoted by \code{col.name} in GRanges object. This
  addition will procede the multiplication. This way the
  bases that have no defined value will have 0 value in the
  Rle vector, and the bases with a predefined value (even
  though it is 0 will be distinguishable in subsequent
  operations)}
}
\value{
  returns a \code{modRleList} object
}
\description{
  This function is similar to \code{coverage} function of
  IRanges and GenomicRanges, however it can deal with
  floating point numbers to some degree. The resulting
  RleList has values from the denoted numeric column, it is
  not simply a coverage vector. In order to deal with
  numbers with decimals some heuristics are applied.
  Warning: Users should NOT use this unless they understand
  what it exactly does.The default arguments will just
  return a coverage vector weighted by the denoted column
  of the GRanges object
}
\note{
  \code{multiply} argument helps users retain some of the
  values with decimal points. Users have to keep in mind
  the fact that \code{add} and/or \code{multiply} arguments
  are in effect in the resulting RleList object and do
  subsequent calculations keep that fact in mind.
}
\seealso{
  \code{\link{coverage}},\code{\link{modRleList}}
}

