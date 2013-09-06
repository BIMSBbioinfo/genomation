genomation: convenience functions for visualing and summarizing genomic intervals
===============
This package is a collection of functions for simplfiying common tasks in genomic feature
analysis. It provides functions for reading BED files as GRanges objects, summarizing genomic
features over predefined windows so users can make average enrichment of features over defined regions
or produce heatmaps.

Reading BED files
______________

read.transcript.features() and read.bed() will read bed files and store them as GRanges
or a list of GRanges objects.



Summarize genomic features(reads, per base scores) over predefined regions 
______________


### Summarize genomic features(reads, per base scores) over predefined regions 
This functionality will summarize given genomic feature over regions of interest.

**target** a RleList or a modRleList or GRanges object to be overlapped with ranges in windows

**windows** a GRanges object that contains the windows of interest. It could be promoters, CpG islands, exons, introns. However the sizes of **windows have to be equal.

**strand.aware** If TRUE (default: FALSE), the strands of the windows will be taken into account in the resulting scoreMatrix. If the strand of a window is -, the values of the bins for that window will be reversed

**ordered** If TRUE (default: FALSE), the input order will be preserved 

*...*  parameters to be passed to modCoverage function. Only needed when target is GRanges.

scoreMatrix(target,windows,bin.num=10,bin.op="mean",strand.aware=FALSE,...)

### Summarize genomic features(reads, per base scores) over bins on predefined regions 
This useful when pre-defined regions are not the same size. This way everything for can be
summarized to bins and then ready to be plotted.


scoreMatrixBin(target,windows,bin.num=10,bin.op="mean",strand.aware=FALSE,...)

Randomize genomic features
______________



Test if overlap between genomic features are due to chance or not
______________

