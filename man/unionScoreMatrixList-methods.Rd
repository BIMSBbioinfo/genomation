\docType{methods}
\name{unionScoreMatrixList}
\alias{unionScoreMatrixList}
\alias{unionScoreMatrixList,ScoreMatrixList-method}
\title{Get common rows from all matrices in a ScoreMatrixList object}
\usage{
  unionScoreMatrixList(sml)

  \S4method{unionScoreMatrixList}{ScoreMatrixList}(sml)
}
\arguments{
  \item{sml}{a \code{ScoreMatrixList} object}
}
\value{
  \code{ScoreMatrixList} object
}
\description{
  Returns a union of rows for each matrix in a
  ScoreMatrixList object. This is done using the rownames
  of each element in the list.
}

