\name{modRleList-class}
\alias{modRleList-class}
\title{An S4 class for storing \code{modCoverage} function results}
\description{
  The resulting object is an extension of a \code{RleList}
  object, and the results of \code{modCoverage} function.
  That is to say, the resulting object is basically an
  \code{RleList} from \code{IRanges} pacakge with extra
  slots defined below.
}
\section{Slots}{
  \describe{ \item{\code{multiply}:}{ number that is
  multiplied with original values used to construct the
  RleList } \item{\code{add}:}{ the value that is added to
  the original values used to construct the RleList} }
}
\seealso{
  \code{\link{modCoverage}}
}

