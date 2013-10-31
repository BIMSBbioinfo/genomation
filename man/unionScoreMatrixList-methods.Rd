\docType{methods}
\name{unionScoreMatrixList}
\alias{unionScoreMatrixList}
\title{Get common rows from a ScoreMatrixList object}
\usage{
  unionScoreMatrixList(sml)
}
\arguments{
  \item{sml}{a \code{ScoreMatrixList} object}
}
\value{
  \code(ScoreMatrixList) object
}
\description{
  Returns a union of rows for each matrix in a
  ScoreMatrixList object. This is done using the rownames
  of each element in the list.
}

