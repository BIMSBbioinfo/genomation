\docType{methods}
\name{randomize.feature}
\alias{randomize.feature}
\alias{randomize.feature,GRanges-method}
\title{function that randomizes the genomic coordinates}
\usage{
  randomize.feature(feature,chrom.sizes=NULL,stranded=TRUE,keep.strand.prop=TRUE,keep.chrom=TRUE,exclude=NULL,include=NULL,seed=NULL)
}
\arguments{
  \item{feature}{a GRanges object to be randomized}

  \item{chrom.sizes}{sizes of chromosomes as a named vector
  (names are chromsomes names and elements of the vectors
  are lengths). , if not given sizes in GRanges object will
  be used if no sizes there the end of each chr will be the
  end last feature on each chr}

  \item{stranded}{if FALSE, all of the returned features
  will be strandless (will have "*" in the strand slot)}

  \item{keep.strand.prop}{If TRUE strands will have the
  same proportion as the features}

  \item{keep.chrom}{If TRUE, number of features and
  randomized features for a chromosome will match.
  Currently seeting this to FALSE is not supported.}

  \item{exclude}{A GRanges object where no randomized
  feature should overlap, can be gaps or unmappable regions
  in the genome as an example.}

  \item{include}{A GRanges object which defines the
  boundaries of randomized features}

  \item{seed}{random number generator seed}
}
\value{
  returns a GRanges object which is randomized version of
  the feature
}
\description{
  This function randomly distributes the coordinates of
  genomic features which is stored in a \code{GRanges}
  object. The randomization can be constrained by supplied
  arguments.
}

