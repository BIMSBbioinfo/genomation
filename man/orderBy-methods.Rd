\docType{methods}
\name{orderBy}
\alias{orderBy}
\alias{orderBy,ScoreMatrixList-method}
\title{Reorder all elements of a ScoreMatrixList to a given ordering vector}
\usage{
  orderBy(sml, ord.vec)
}
\arguments{
  \item{sml}{\code{ScoreMatrixList} object}

  \item{ord.vec}{an integer vector}
}
\value{
  \code{ScoreMatrixList} object
}
\description{
  Reorder all elements of a ScoreMatrixList to a given
  ordering vector
}
\examples{
data(cage)
 data(cpgi)
 data(promoters)

 cage$tpm = NULL
 target = GRangesList(cage=cage, cpgi=cpgi)
 sml = ScoreMatrixList(target, promoters, bin.num=10)
 kmeans.clust = kmeans(sml$cage,3)

 sml.ordered = orderBy(sml, kmeans.clust$cluster)
 multiHeatMatrix(sml.ordered)
}

