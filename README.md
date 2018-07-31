<a name="logo"/>
<div align="center">
<img src="https://github.com/BIMSBbioinfo/genomation/blob/master/inst/genomation_logo.png" alt="hex Logo"  ></img>
</a>
</div>

# genomation: convenience functions for visualing and summarizing genomic intervals

Status [![Build Status](https://api.travis-ci.org/BIMSBbioinfo/genomation.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/genomation) [![codecov.io](https://codecov.io/github/BIMSBbioinfo/genomation/coverage.svg?branch=master)](https://codecov.io/github/BIMSBbioinfo/genomation?branch=master)    [![BioC_years](http://www.bioconductor.org/shields/years-in-bioc/genomation.svg)](http://www.bioconductor.org/packages/release/bioc/html/genomation.html) [![BioC_availability](http://www.bioconductor.org/shields/availability/release/genomation.svg)](http://www.bioconductor.org/packages/release/bioc/html/genomation.html)


This package is a collection of functions for simplfiying common tasks in genomic feature/interval
analysis. It provides functions for reading BED and GFF files as GRanges objects, summarizing genomic features over predefined windows so users can make average enrichment of features over defined regions or produce heatmaps. It can also annotate given regions
with other genomic features such as exons,introns and promoters.

# Installation

### install via Bioconductor
```R
if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
BiocManager::install("genomation")

```

### Install the latest version via install_github
```R
#' Install dependecies
install.packages( c("data.table","plyr","reshape2","ggplot2","gridBase","devtools"))
if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges","rtracklayer","impute","Rsamtools"))


#' install the packages
library(devtools)
install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)

```

# Using the package
see the package vignette [here](http://bioconductor.org/packages/release/bioc/vignettes/genomation/inst/doc/GenomationManual.html)
