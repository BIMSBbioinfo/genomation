\docType{methods}
\name{annotateWithGeneParts}
\alias{annotateWithGeneParts}
\alias{annotateWithGeneParts,GRanges,GRangesList-method}
\alias{annotateWithGeneParts,GRangesList,GRangesList-method}
\title{Annotate given  object with promoter, exon, intron and intergenic regions}
\usage{
  annotateWithGeneParts(target, GRangesList.obj,
    strand = F)
}
\arguments{
  \item{target}{A \code{\link{GRanges}} or
  \code{\link{GRangesList}} object storing chromosome
  locations to be annotated (e.g. chipseq peaks)}

  \item{GRangesList.obj}{A GRangesList object containing
  GRanges object for promoter, exons, introns and TSSes, or
  simply output of readTranscriptFeatures function}

  \item{strand}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}
}
\value{
  \code{annotationByGenicParts} object or a list of
  \code{annotationByGenicParts} objects if target is a
  \code{\link{GRangesList}} object.
}
\description{
  The function annotates \code{GRangesList} or
  \code{GRanges} object as overlapping with
  promoter,exon,intron or intergenic regions.
}

