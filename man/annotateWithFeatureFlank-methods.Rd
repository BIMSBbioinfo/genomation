\docType{methods}
\name{annotateWithFeatureFlank}
\alias{annotateWithFeatureFlank}
\alias{annotateWithFeatureFlank,GRanges,GRanges,GRanges-method}
\title{Function to annotate a given GRanges object with promoter,exon,intron & intergenic values}
\usage{
annotateWithFeatureFlank(target,feature,flank,feature.name="feat",
                                flank.name="flank",strand=FALSE)

\S4method{annotateWithFeatureFlank}{GRanges,GRanges,GRanges}(target, feature,
  flank, feature.name = NULL, flank.name = "flank", strand = FALSE,
  intersect.chr = FALSE)
}
\arguments{
  \item{target}{a granges object storing chromosome
  locations to be annotated}

  \item{feature}{a granges object storing chromosome
  locations of a feature (can be CpG islands, ChIP-seq
  peaks, etc)}

  \item{flank}{a granges object storing chromosome
  locations of the flanks of the feature}

  \item{feature.name}{string for the name of the feature}

  \item{flank.name}{string for the name of the flanks}

  \item{strand}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FAULT)}

  \item{intersect.chr}{boolean, whether to select only
  chromosomes that are common to feature and target. FALSE
  by default}
}
\value{
  returns an \code{AnnotationByFeature} object
}
\description{
  Function to annotate a given GRanges object with
  promoter,exon,intron & intergenic values
}
\examples{
data(cpgi)
data(cage)
cpgi.flanks = getFlanks(cpgi)
flank.annot = annotateWithFeatureFlank(cage, cpgi, cpgi.flanks)
}

