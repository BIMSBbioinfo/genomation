\docType{methods}
\name{readBroadPeak}
\alias{readBroadPeak}
\title{A function to read the Encode formatted broad peak file into a GRanges object}
\usage{
  readBroadPeak(file)
}
\arguments{
  \item{file}{a abosulte or relative path to a bed file
  formatted by the Encode broadPeak standard}
}
\value{
  a GRanges object
}
\description{
  A function to read the Encode formatted broad peak file
  into a GRanges object
}
\examples{
library(genomationData)
broad.peak.file = list.files(system.file('extdata', package='genomationData'),
                             full.names=T, pattern='broadPeak')[1]
broad.peak = readBroadPeak(broad.peak.file)
head(broad.peak)
}

