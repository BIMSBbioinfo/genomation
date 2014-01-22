\docType{methods}
\name{randomizeFeature}
\alias{randomizeFeature}
\alias{randomizeFeature,GRanges-method}
\title{function that randomizes the genomic coordinates}
\usage{
randomizeFeature(feature, chrom.sizes = NULL, stranded = TRUE,
  keep.strand.prop = TRUE, keep.chrom = TRUE, exclude = NULL,
  include = NULL, seed = NULL, nrand = 1)

\S4method{randomizeFeature}{GRanges}(feature, chrom.sizes, stranded,
  keep.strand.prop, keep.chrom, exclude, include, seed, nrand)
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
  boundaries of randomized features. If not provided the
  whole genome is used, as defined using the chrom.sizes
  parameter.}

  \item{seed}{random number generator seed}

  \item{nrand}{number of randomizations (default:1)}
}
\value{
returns a GRanges object which is randomized version of the
feature, along with a "set" column in the metadata which
designates to which iteration of the randomization the
range belong.
}
\description{
This function randomly distributes the coordinates of
genomic features which is stored in a \code{GRanges}
object. The randomization can be constrained by supplied
arguments. The function is still in Beta mode - the regions
can overlap excluded regions, and the randomized regions
are not disjoint. Please take care that the excluded and
included regions are not too strict when compared to the
total width of the ranges.
}

