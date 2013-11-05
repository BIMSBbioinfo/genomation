\docType{methods}
\name{annotateWithFeature}
\alias{annotateWithFeature}
\alias{annotateWithFeature,GRanges,GRanges-method}
\title{Function to annotate given GRanges object with a given genomic feature}
\usage{
  annotateWithFeature(target,feature,strand=FALSE,extend=0,feature.name="feat1")

  \S4method{annotateWithFeature}{GRanges,GRanges}(target,
    feature, strand = FALSE, extend = 0,
    feature.name = "feat1")
}
\arguments{
  \item{target}{a GRanges object storing chromosome
  locations to be annotated}

  \item{feature}{a GRanges object storing chromosome
  locations of a feature (can be CpG islands, ChIP-seq
  peaks, etc)}

  \item{strand}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FAULT)}

  \item{extend}{specifiying a positive value will extend
  the feature on both sides as much as \code{extend}}

  \item{feature.name}{name of the annotation feature. For
  example: H3K4me1,CpGisland etc.}
}
\value{
  returns an \code{annotationByFeature} object
}
\description{
  Function to annotate given GRanges object with a given
  genomic feature
}

