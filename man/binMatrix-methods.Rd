\docType{methods}
\name{binMatrix}
\alias{binMatrix}
\title{Bins the columns of a matrix using a user provided function}
\usage{
  binMatrix(mat, nbins = NULL, fun = "mean", ...)
}
\arguments{
  \item{mat}{a \code{scoreMatrix} object}

  \item{nbins}{a \code{integer} number of bins in the final
  matrix}

  \item{fun}{a \code{character} vector representing the
  function to be used for bining}
}
\value{
  \code{scoreMatrix} object
}
\description{
  Bins the columns of a matrix using a user provided
  function
}

