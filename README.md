# genomation: convenience functions for visualing and summarizing genomic intervals

This package is a collection of functions for simplfiying common tasks in genomic feature
analysis. It provides functions for reading BED and GFF files as GRanges objects, summarizing genomic features over predefined windows so users can make average enrichment of features over defined regions or produce heatmaps. It can also annotate given regions
with other genomic features such as exons,introns and promoters.

# Installation

```R
#' Install dependecies
install.packages( c("data.table","plyr","reshape2","ggplot2","gridBase","devtools"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","rtracklayer","impute","Rsamtools"))

#' install the packages
library(devtools)
install_github("genomation", username = "al2na")

}
```

# Using the package
see the package vignette [here](https://github.com/al2na/genomation/raw/development/inst/doc/GenomationManual-knitr.pdf)
