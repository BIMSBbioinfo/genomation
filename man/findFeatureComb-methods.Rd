\docType{methods}
\name{findFeatureComb}
\alias{findFeatureComb}
\alias{findFeatureComb,GRanges,GRangesList,GRangesList-methods}
\title{Find combitations of genomic features}
\arguments{
  \item{gl}{a \code{GRangesList} object, containing ranges
  for which represent regions enriched for transcription
  factor binding}

  \item{width}{\code{integer} is the requested width of
  each enriched region. If 0 the ranges are not resized, if
  a positive integer, the width of all ranges is set to
  that number. Ranges are resized relative to the center of
  original ranges.}

  \item{use.names}{a boolean which tells the function
  whether to return the resulting ranges with a numeric
  vector which designates each class (the default), or to
  construct the names of each class using the names from
  the GRangesList}

  \item{collapse.char}{a character which will be used to
  separate the class names if use.names=TRUE. The default
  is ':'}
}
\value{
  a GRanges object
}
\description{
  Provided a GRangesList, finds the combinations of sets of
  ranges. It is mostly used to look at the combinatorics of
  transcription factor binding. The function works by,
  firstly, constructing a union of all ranges in the list,
  which are then designated by the combinatorics of overlap
  with the original sets. A caveat of this approach is that
  the number of possible combinations increases
  exponentially, so we would advise you to use it with up
  to 6 data sets. If you wish to take a look at a greater
  number of factors, methods like self organizing maps or
  ChromHMM might be more appropriate.
}

