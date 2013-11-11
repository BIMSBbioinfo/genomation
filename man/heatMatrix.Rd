\name{heatMatrix}
\alias{heatMatrix}
\title{Draw a heatmap of a given ScoreMatrix object}
\usage{
  heatMatrix(mat, grid = FALSE, col = NULL, xcoords = NULL,
    group = NULL, group.col = NULL, order = FALSE,
    winsorize = c(0, 100), kmeans = FALSE, k = 3,
    main = "", legend.name = NULL, cex.legend = 1,
    xlab = NULL, cex.main = 1, cex.lab = 1, cex.axis = 1,
    newpage = TRUE)
}
\arguments{
  \item{mat}{a \code{ScoreMatrix} object}

  \item{grid}{if TRUE, grid graphics will be used. if
  FALSE, base graphics will be used on the top level, so
  users can use par(mfrow) or par(mfcol) prior to calling
  the function. Default:FALSE}

  \item{col}{a vector of colors, such as the ones created
  by heat.colors(10).  If NULL (which is default), jet
  color scheme (common in matlab plots) will be used.}

  \item{xcoords}{a vector of numbers showing relative
  positions of the bases or windows. It must match the
  number of columns in the \code{ScoreMatrix}.
  Alternatively, it could be a numeric vector of two
  elements. Such as c(0,100) showing the relative start and
  end coordinates of the first and last column of the
  \code{ScoreMatrix} object.}

  \item{group}{a list of vectors of row numbers or a
  factor. This grouping is used for rowside colors of the
  heatmap. If it is a list, each element of the list must
  be a vector of row numbers. Names of the elements of the
  list will be used as names of groups.  If \code{group} is
  a factor , it's length must match the number of rows of
  the matrix, and factor levels will be used as the names
  of the groups in the plot.}

  \item{group.col}{a vector of color names to be used at
  the rowside colors if \code{group} argument is given or
  \code{kmeans=TRUE}}

  \item{order}{Logical indicating if the rows should be
  ordered or not (Default:FALSE). If \code{order=TRUE} the
  matrix will be ordered with rowSums(mat) values in
  descending order. If kmeans=TRUE or \code{group} argument
  is provided, first the groups/clusters will be ordered in
  descending order of sums of rows then, everything within
  the clusters will be ordered by sums of rows.}

  \item{winsorize}{Numeric vector of two, defaults to
  c(0,100). This vector determines the upper and lower
  percentile values to limit the extreme values. For
  example, c(0,99) will limit the values to only 99th
  percentile, everything above the 99 percentile will be
  equalized to the value of 99th percentile.This is useful
  for visualization of matrices that have outliers.}

  \item{kmeans}{Logical indicating if kmeans clustering
  should be done on the rows or not (Default:FALSE).}

  \item{k}{Defaults to 3. It designates the number of
  clusters to be returned by kmeans clustering.}

  \item{main}{a character string for the plot title}

  \item{legend.name}{a character label plotted next to the
  legend}

  \item{cex.legend}{A numerical value giving the amount by
  which legend axis marks should be magnified relative to
  the default}

  \item{xlab}{label a character string for x-axis of the
  heatmap}

  \item{cex.lab}{A numerical value giving the amount by
  which axis labels (including 'legend.name') should be
  magnified relative to the default.}

  \item{cex.main}{A numerical value giving the amount by
  which plot title should be magnified}

  \item{cex.axis}{A numerical value giving the amount by
  which axis marks should be magnified relative to the
  default}

  \item{newpage}{logical indicating if
  \code{grid.newpage()} function should be invoked if
  \code{grid=TRUE}.}
}
\value{
  returns kmeans clustering result invisibly, if
  kmeans=TRUE
}
\description{
  The function makes a heatmap out of given
  \code{ScoreMatrix} object. If desired it can use
  clustering using k-means and plot cluster color codes as
  a sidebar. In addition, user can define groups of rows
  using 'group' argument.
}
\examples{
data(cage)
data(promoters)
scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
                    weight.col="tpm")



heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
          xcoords=-1000:1000,
          cex.legend=0.8,main="GAGE clusters on promoters",cex.lab=1,
          cex.axis=0.9,grid=FALSE)

set.seed(1000)
heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
          xcoords=-1000:1000,kmeans=TRUE,k=3,
          cex.legend=0.8,main="GAGE clusters on promoters",cex.lab=1,
          cex.axis=0.9,grid=FALSE)
}

