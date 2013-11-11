\docType{methods}
\name{plotGenicAnnotation}
\alias{plotGenicAnnotation}
\alias{plotGenicAnnotation,annotationByFeature-method}
\title{Plots the enrichment of each feature in the set in the gene annotation}
\usage{
  plotGenicAnnotation(l, cluster = FALSE,
    col = c("white", "cornflowerblue"))
}
\arguments{
  \item{l}{a \code{list} of \code{annotationByGenicParts}
  objects}

  \item{cluster}{TRUE/FALSE. If TRUE the heatmap is going
  to be clustered using hierarchical clustering}

  \item{col}{a vector of two colors that will be used for
  interpolation.  The first color is the lowest one, the
  second is the highest one usage plotGenicAnnotation(l,
  cluster=FALSE, col=c('white','cornflowerblue'),...)}
}
\description{
  This function plots a heatmap of enrichment of each range
  in given gene feature
}

