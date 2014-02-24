\docType{methods}
\name{getAssociationWithTSS}
\alias{AnnotationByGeneParts-method}
\alias{getAssociationWithTSS}
\alias{getAssociationWithTSS,}
\alias{getAssociationWithTSS,-methods}
\alias{getAssociationWithTSS,AnnotationByGeneParts-method}
\title{Get distance to nearest TSS and gene id from AnnotationByGeneParts}
\usage{
getAssociationWithTSS(x)

\S4method{getAssociationWithTSS}{AnnotationByGeneParts}(x)
}
\arguments{
  \item{x}{a \code{AnnotationByGeneParts} object}
}
\value{
RETURNS a data.frame containing row number of the target
features, distance of target to nearest TSS, TSS/Gene name,
TSS strand
}
\description{
This accessor function gets the nearest TSS, its distance
to target feature, strand and name of TSS/gene from
AnnotationByGeneParts object
}
\examples{
data(cage)
bed.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
head(getAssociationWithTSS(cage.annot))
}

