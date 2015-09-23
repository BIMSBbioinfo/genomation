(require("genomation") & require("GenomicRanges") & require("GenomicAlignments") & require("rtracklayer") & require("readr") & require("Biostrings") & require("BSgenome")) || stop("unable to load genomation package")
genomation:::.test()
