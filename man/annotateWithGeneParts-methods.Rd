\docType{methods}
\name{annotateWithGeneParts}
\alias{annotateWithGeneParts}
\alias{annotateWithGeneParts,GRanges,GRangesList-method}
\alias{annotateWithGeneParts,GRangesList,GRangesList-method}
\title{Annotate given  object with promoter, exon, intron and intergenic regions}
\usage{
  annotateWithGeneParts(target, feature, strand = F,
    intersect.chr = FALSE)

  \S4method{annotateWithGeneParts}{GRanges,GRangesList}(target,
    feature, strand = F, intersect.chr = FALSE)

  \S4method{annotateWithGeneParts}{GRangesList,GRangesList}(target,
    feature, strand = F, intersect.chr = FALSE)
}
\arguments{
  \item{target}{\code{\link{GRanges}} or
  \code{\link{GRangesList}} object storing chromosome
  locations to be annotated (e.g. chipseq peaks)}

  \item{feature}{GRangesList object containing GRanges
  object for promoter, exons, introns and transcription
  start sites, or simply output of readTranscriptFeatures
  function}

  \item{strand}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}

  \item{intersect.chr}{boolean, whether to select only
  chromosomes that are common to feature and target. FALSE
  by default}
}
\value{
  \code{AnnotationByGeneParts} object or a list of
  \code{AnnotationByGeneParts} objects if target is a
  \code{\link{GRangesList}} object.
}
\description{
  The function annotates \code{GRangesList} or
  \code{GRanges} object as overlapping with
  promoter,exon,intron or intergenic regions.
}
\examples{
data(cage)
bed.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
}

