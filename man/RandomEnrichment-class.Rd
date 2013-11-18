\docType{class}
\name{RandomEnrichment-class}
\alias{RandomEnrichment-class}
\title{An S4 class for storing \code{getRandomEnrichment} function results}
\description{
  The resulting object stores the results of
  \code{getRandomEnrichment} function
}
\section{Slots}{
  \describe{ \item{\code{orig.cnt}:}{number of features
  overlapping with query at \code{getRandomEnrichment} }
  \item{\code{rand.olap.dist}:}{set of number of features
  overlapping with randomized queries at
  \code{getRandomEnrichment}} \item{\code{log2fc}:}{ log2
  fold change calculated by dividing \code{orig.cnt} by
  mean(\code{rand.olap.dist}) and taking log2 of that
  result} \item{\code{p.value}:}{P-value assuming
  \code{rand.olap.dist} has a normal distribution and
  comparing \code{orig.cnt} with that distribution }
  \item{\code{rand.p.value}:}{ p-value from randomization
  by calculation the proportion of how many times a random
  number of overlap exceeds the original number of overlap}
  }
}
\seealso{
  \code{\link{getRandomEnrichment}}
}

