\docType{methods}
\name{annotateWithGeneParts}
\alias{annotateWithGeneParts}
\alias{annotateWithGeneParts,GRanges,GRangesList-method}
\title{Function to annotate given GRanges object with promoter, exon, intron & intergenic ranges}
\arguments{
  \item{target:}{A GRanges object storing chromosome
  locations to be annotated (e.g. chipseq peaks)}

  \item{GRangesList.obj:}{A GRangesList object containing
  GRanges object for promoter, exons, introns and TSSes, or
  simply output of readTranscriptFeatures function}

  \item{strand:}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}

  \item{target:}{A GRangesList object storing chromosome
  locations to be annotated (e.g. chipseq peaks from
  multiple experiments)}

  \item{GRangesList.obj:}{A GRangesList object containing
  GRanges object for promoter, exons, introns and TSSes, or
  simply output of readTranscriptFeatures function}

  \item{strand:}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}
}
\value{
  \code{annotationByGenicParts} object

  \code{annotationByGenicParts} object
}
\description{
  Function to annotate given GRanges object with promoter,
  exon, intron & intergenic ranges

  Given a GRangesList object it annotates each Range with
  gene annotation
}

