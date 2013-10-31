\docType{methods}
\name{heatMatrix}
\alias{heatMatrix}
\alias{heatMatrix,scoreMatrix-method}
\title{visual representation of scoreMatrix using a heatmap
The rows can be reordered using one factor and one numeric vector}
\usage{
  heatMatrix(mat, fact = NULL, add.sep = TRUE,
    ord.vec = NULL, shift = 0, mat.cols = NULL,
    fact.cols = NULL, xlab = "Position", ylab = "Region",
    main = "Positional profile", class.names = NULL,
    use.names = FALSE, ...)
}
\arguments{
  \item{mat}{a \code{scoreMatrix} object}

  \item{fact}{a \code{factor} of length equal to
  \code{nrow(mat)}. Unused factor levels are dropped}

  \item{ord.vec}{a \code{vector} of class \code{numeric} of
  the same length as mat, which is going to be used for
  ordering of the rows}

  \item{shift}{shift the start coordinate of the x axis
  (plot starts at -shift)}

  \item{mat.cols}{a vector of colors used for plotting of
  the heatmap. Default colors range from lightgray to
  darkblue.}

  \item{fact.cols}{a vector of colors used for plotting of
  the factor key}

  \item{xlab}{x axis label}

  \item{ylab}{y axis label}

  \item{main}{plot name}

  \item{class.names}{names for each factor class - has to
  have the same lenght as \code{levels(fact)}}

  \item{...}{other options to be passed to functions
  (obsolete at the moment)}
}
\value{
  nothing
}
\description{
  visual representation of scoreMatrix using a heatmap The
  rows can be reordered using one factor and one numeric
  vector
}
\examples{
data(cage)
  data(promoters)
  myMat2=scoreMatrix(target=cage,windows=promoters,
                        col.name="tpm",strand.aware=TRUE)
  plot(colMeans(myMat2,na.rm=TRUE),type="l")
  heatMatrix(myMat2,fact=)
}

