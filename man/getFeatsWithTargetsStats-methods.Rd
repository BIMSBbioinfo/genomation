\docType{methods}
\name{getFeatsWithTargetsStats}
\alias{getFeatsWithTargetsStats}
\alias{getFeatsWithTargetsStats,AnnotationByFeature-method}
\title{Get the percentage/count of annotation features overlapping with target
features from AnnotationByFeature}
\usage{
getFeatsWithTargetsStats(x,percentage=TRUE)

\S4method{getFeatsWithTargetsStats}{AnnotationByFeature}(x, percentage)
}
\arguments{
  \item{x}{a \code{AnnotationByFeature} object}

  \item{percentage}{TRUE|FALSE. If TRUE percentage of
  annotation features will be returned. If FALSE, number of
  annotation features will be returned}
}
\value{
RETURNS a vector of percentages or counts showing quantity
of annotation features overlapping with target features
}
\description{
This function retrieves percentage/number of annotation
features overlapping with targets. For example, if
\code{AnnotationByFeature} object is containing statistics
of differentially methylated regions overlapping with gene
annotation. This function will return number/percentage of
introns,exons and promoters overlapping with differentially
methylated regions.
}
\examples{
data(cage)
bed.file=system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
getFeatsWithTargetsStats(cage.annot)
}

