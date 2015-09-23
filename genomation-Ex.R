pkgname <- "genomation"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('genomation')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AnnotationByGeneParts-methods")
### * AnnotationByGeneParts-methods

flush(stderr()); flush(stdout())

### Name: getAssociationWithTSS
### Title: Get distance to nearest TSS and gene id from
###   AnnotationByGeneParts
### Aliases: AnnotationByGeneParts-method getAssociationWithTSS
###   getAssociationWithTSS, getAssociationWithTSS,-methods
###   getAssociationWithTSS,AnnotationByGeneParts-method

### ** Examples

data(cage)
bed.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
head(getAssociationWithTSS(cage.annot))



cleanEx()
nameEx("ScoreMatrix-methods")
### * ScoreMatrix-methods

flush(stderr()); flush(stdout())

### Name: ScoreMatrix
### Title: Get base-pair score for bases in each window
### Aliases: ScoreMatrix ScoreMatrix,GRanges,GRanges-method
###   ScoreMatrix,RleList,GRanges-method
###   ScoreMatrix,character,GRanges-method

### ** Examples

# When target is GRanges
data(cage)
data(promoters)
scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
                                weight.col="tpm")


# When target is RleList
library(GenomicRanges)
covs = coverage(cage)
scores2 = ScoreMatrix(target=covs,windows=promoters,strand.aware=TRUE)
scores2

# When target is a bam file
bam.file = system.file('unitTests/test.bam', package='genomation')
windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5))
scores3 = ScoreMatrix(target=bam.file,windows=windows, type='bam')
scores3



cleanEx()
nameEx("ScoreMatrixBin-methods")
### * ScoreMatrixBin-methods

flush(stderr()); flush(stdout())

### Name: ScoreMatrixBin
### Title: Get bin score for bins on each window
### Aliases: ScoreMatrixBin ScoreMatrixBin,GRanges,GRanges-method
###   ScoreMatrixBin,RleList,GRanges-method
###   ScoreMatrixBin,character,GRanges-method

### ** Examples

data(cage)
data(cpgi)
data(promoters)
myMat=ScoreMatrixBin(target=cage,
                      windows=cpgi,bin.num=10,bin.op="mean",weight.col="tpm")
## No test: 
plot(colMeans(myMat,na.rm=TRUE),type="l")
## End(No test)

myMat2=ScoreMatrixBin(target=cage,
                       windows=promoters,bin.num=10,bin.op="mean",
                       weight.col="tpm",strand.aware=TRUE)
## No test: 
plot(colMeans(myMat2,na.rm=TRUE),type="l")
## End(No test)



cleanEx()
nameEx("ScoreMatrixList-methods")
### * ScoreMatrixList-methods

flush(stderr()); flush(stdout())

### Name: ScoreMatrixList
### Title: Make ScoreMatrixList from multiple targets
### Aliases: ScoreMatrixList

### ** Examples

# visualize the distribution of cage clusters and cpg islands around promoters
library(GenomicRanges)
data(cage)
data(cpgi)
data(promoters)

cage$tpm = NULL
targets = GRangesList(cage=cage, cpgi=cpgi)
sml = ScoreMatrixList(targets, promoters, bin.num=10, strand.aware=TRUE)
sml
## No test: 
multiHeatMatrix(sml)
## End(No test)



cleanEx()
nameEx("annotateWithFeature-methods")
### * annotateWithFeature-methods

flush(stderr()); flush(stdout())

### Name: annotateWithFeature
### Title: Function to annotate given GRanges object with a given genomic
###   feature
### Aliases: annotateWithFeature annotateWithFeature,GRanges,GRanges-method

### ** Examples

data(cpgi)
data(promoters)
annot = annotateWithFeature(cpgi, promoters)



cleanEx()
nameEx("annotateWithFeatureFlank-methods")
### * annotateWithFeatureFlank-methods

flush(stderr()); flush(stdout())

### Name: annotateWithFeatureFlank
### Title: Function to annotate a given GRanges object with
###   promoter,exon,intron & intergenic values
### Aliases: annotateWithFeatureFlank
###   annotateWithFeatureFlank,GRanges,GRanges,GRanges-method

### ** Examples

data(cpgi)
data(cage)
cpgi.flanks = getFlanks(cpgi)
flank.annot = annotateWithFeatureFlank(cage, cpgi, cpgi.flanks)



cleanEx()
nameEx("annotateWithGeneParts-methods")
### * annotateWithGeneParts-methods

flush(stderr()); flush(stdout())

### Name: annotateWithGeneParts
### Title: Annotate given object with promoter, exon, intron and intergenic
###   regions
### Aliases: annotateWithGeneParts
###   annotateWithGeneParts,GRanges,GRangesList-method
###   annotateWithGeneParts,GRangesList,GRangesList-method

### ** Examples

data(cage)
bed.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
cage.annot



cleanEx()
nameEx("binMatrix-methods")
### * binMatrix-methods

flush(stderr()); flush(stdout())

### Name: binMatrix
### Title: Bins the columns of a matrix using a user provided function
### Aliases: binMatrix binMatrix,ScoreMatrix-method
###   binMatrix,ScoreMatrixList-method

### ** Examples

# binning the columns in a ScoreMatrix object
library(GenomicRanges)
target = GRanges(rep(c(1,2),each=7), IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
weight = rep(c(1,2),each=7),
strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))
windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5),
strand=c('-','+','-','+'))
sm = ScoreMatrix(target, windows)
bin = binMatrix(sm, bin.num=2)



cleanEx()
nameEx("convertBed2Exons-methods")
### * convertBed2Exons-methods

flush(stderr()); flush(stdout())

### Name: convertBed2Exons
### Title: convert a data frame read-in from a bed file to a GRanges object
###   for exons
### Aliases: convertBed2Exons convertBed2Exons,data.frame-method

### ** Examples

file = system.file('extdata/chr21.refseq.hg19.bed', package='genomation')
bed12 = read.table(file)
exons = convertBed2Exons(bed12)
head(exons)



cleanEx()
nameEx("convertBed2Introns-methods")
### * convertBed2Introns-methods

flush(stderr()); flush(stdout())

### Name: convertBed2Introns
### Title: convert a data frame read-in from a bed file to a GRanges object
###   for introns
### Aliases: convertBed2Introns convertBed2Introns,data.frame-method

### ** Examples

file = system.file('extdata/chr21.refseq.hg19.bed', package='genomation')
bed12 = read.table(file)
introns = convertBed2Introns(bed12)
head(introns)



cleanEx()
nameEx("findFeatureComb-methods")
### * findFeatureComb-methods

flush(stderr()); flush(stdout())

### Name: findFeatureComb
### Title: Find combitations of genomic features
### Aliases: findFeatureComb findFeatureComb,GRangesList-method

### ** Examples

library(GenomicRanges)
g = GRanges(paste('chr',rep(1:2, each=3), sep=''), IRanges(rep(c(1,5,9), times=2), width=3))
gl = GRangesList(g1=g, g2=g[2:5], g3=g[3:4])
findFeatureComb(gl)
findFeatureComb(gl, use.names=TRUE)



cleanEx()
nameEx("getFeatsWithTargetsStats-methods")
### * getFeatsWithTargetsStats-methods

flush(stderr()); flush(stdout())

### Name: getFeatsWithTargetsStats
### Title: Get the percentage/count of annotation features overlapping with
###   target features from AnnotationByFeature
### Aliases: getFeatsWithTargetsStats
###   getFeatsWithTargetsStats,AnnotationByFeature-method

### ** Examples

data(cage)
bed.file=system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
getFeatsWithTargetsStats(cage.annot)



cleanEx()
nameEx("getFlanks-methods")
### * getFlanks-methods

flush(stderr()); flush(stdout())

### Name: getFlanks
### Title: Function to get upstream and downstream adjecent regions to a
###   genomic feature such as CpG islands
### Aliases: getFlanks getFlanks,GRanges-method

### ** Examples

data(cpgi)
cpgi.flanks = getFlanks(cpgi)
head(cpgi.flanks)



cleanEx()
nameEx("getRandomEnrichment-methods")
### * getRandomEnrichment-methods

flush(stderr()); flush(stdout())

### Name: getRandomEnrichment
### Title: get enrichment based on randomized feature overlap
### Aliases: getRandomEnrichment getRandomEnrichment,GRanges,GRanges-method

### ** Examples

data(cage)
data(cpgi)
## No test: 
enr = getRandomEnrichment(cage, cpgi, randomizations=50)
## End(No test)



cleanEx()
nameEx("getTargetAnnotationStats-methods")
### * getTargetAnnotationStats-methods

flush(stderr()); flush(stdout())

### Name: getTargetAnnotationStats
### Title: Get the percentage of target features overlapping with
###   annotation from AnnotationByFeature
### Aliases: getTargetAnnotationStats
###   getTargetAnnotationStats,AnnotationByFeature-method

### ** Examples

data(cage)
bed.file=system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
cage.annot=annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
getTargetAnnotationStats(cage.annot)



cleanEx()
nameEx("gffToGRanges")
### * gffToGRanges

flush(stderr()); flush(stdout())

### Name: gffToGRanges
### Title: Converts a gff formated data.frame into a GenomicRanges object.
###   The GenomicRanges object needs to be properly formated for the
###   function to work.
### Aliases: gffToGRanges

### ** Examples

gff.file = system.file('extdata/chr21.refseq.hg19.gtf', package='genomation')
gff = gffToGRanges(gff.file)



cleanEx()
nameEx("heatMatrix")
### * heatMatrix

flush(stderr()); flush(stdout())

### Name: heatMatrix
### Title: Draw a heatmap of a given ScoreMatrix object
### Aliases: heatMatrix

### ** Examples

data(cage)
data(promoters)
scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
                   weight.col="tpm")

set.seed(1000)
## No test: 
heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
           xcoords=-1000:1000,
           cex.legend=0.8,main="CAGE clusters on promoters",cex.lab=1,
           cex.axis=0.9,grid=FALSE)

## examples using clustering functions
## k-means
cl1 <- function(x) kmeans(x, centers=3)$cluster
set.seed(1000)
heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
         xcoords=-1000:1000,clustfun=cl1,
         cex.legend=0.8,main="CAGE clusters on promoters",cex.lab=1,
         cex.axis=0.9,grid=FALSE,
         user.order=c(1,3,2))

## hierarchical clustering
cl2 <- function(x) cutree(hclust(dist(x), method="complete"), k=3)
set.seed(1000)
heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
         xcoords=-1000:1000,clustfun=cl2,
         cex.legend=0.8,main="CAGE clusters on promoters",cex.lab=1,
         cex.axis=0.9,grid=FALSE)
## End(No test)



cleanEx()
nameEx("heatMeta")
### * heatMeta

flush(stderr()); flush(stdout())

### Name: heatMeta
### Title: Heatmap for meta-region profiles
### Aliases: heatMeta

### ** Examples

data(cage)
 data(promoters)
 scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE)
 data(cpgi)
 scores2=ScoreMatrix(target=cpgi,windows=promoters,strand.aware=TRUE)

 x=new("ScoreMatrixList",list(scores1,scores2))
 ## No test: 
 heatMeta(mat=x,legend.name="fg",cex.legend=0.8,main="fdf",cex.lab=6,
          cex.axis=0.9)
 
## End(No test)



cleanEx()
nameEx("intersectScoreMatrixList-methods")
### * intersectScoreMatrixList-methods

flush(stderr()); flush(stdout())

### Name: intersectScoreMatrixList
### Title: Get common rows from all matrices in a ScoreMatrixList object
### Aliases: intersectScoreMatrixList
###   intersectScoreMatrixList,ScoreMatrixList-method

### ** Examples

library(GenomicRanges)
target = GRanges(rep(c(1,2),each=7),
                  IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
                  weight = rep(c(1,2),each=7))

windows1 = GRanges(rep(c(1,2),each=2),
                    IRanges(rep(c(1,2), times=2), width=5),
                    strand=c('-','+','-','+'))
windows2 = windows1[c(1,3)]
sml = as(list(ScoreMatrix(target, windows1),
               ScoreMatrix(target, windows2)), 'ScoreMatrixList')
sml
intersectScoreMatrixList(sml)



cleanEx()
nameEx("multiHeatMatrix")
### * multiHeatMatrix

flush(stderr()); flush(stdout())

### Name: multiHeatMatrix
### Title: Draw multiple heatmaps from a ScoreMatrixList object
### Aliases: multiHeatMatrix

### ** Examples

data(cage)
data(promoters)
scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE)

data(cpgi)
scores2=ScoreMatrix(target=cpgi,windows=promoters,strand.aware=TRUE)

sml=new("ScoreMatrixList",list(a=scores1,b=scores2))

# use with k-means
## No test: 
multiHeatMatrix(sml,
                 clustfun=function(x) kmeans(x, centers=2)$cluster,
                 cex.axis=0.8,xcoords=c(-1000,1000),
                 winsorize=c(0,99),
                 legend.name=c("tpm","coverage"),xlab="region around TSS")

# use with hierarchical clustering
cl2 <- function(x) cutree(hclust(dist(x), method="complete"), k=2)
multiHeatMatrix(sml,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
         xcoords=-1000:1000,clustfun=cl2,
         cex.legend=0.8,cex.lab=1,
         cex.axis=0.9,grid=FALSE)

# use different colors
require(RColorBrewer)
col.cage= brewer.pal(9,"Blues")
col.cpgi= brewer.pal(9,"YlGn")
multiHeatMatrix(sml,
                 clustfun=function(x) kmeans(x, centers=2)$cluster,
                 cex.axis=0.8,xcoords=c(-1000,1000),
                 winsorize=c(0,99),col=list(col.cage,col.cpgi),
                 legend.name=c("tpm","coverage"),xlab="region around TSS")
## End(No test)



cleanEx()
nameEx("orderBy-methods")
### * orderBy-methods

flush(stderr()); flush(stdout())

### Name: orderBy
### Title: Reorder all elements of a ScoreMatrixList to a given ordering
###   vector
### Aliases: orderBy orderBy,ScoreMatrixList-method

### ** Examples

library(GenomicRanges)
data(cage)
data(cpgi)
data(promoters)

cage$tpm = NULL
targets = GRangesList(cage=cage, cpgi=cpgi)
sml = ScoreMatrixList(targets, promoters, bin.num=10)
kmeans.clust = kmeans(sml$cage,3)

sml.ordered = orderBy(sml, kmeans.clust$cluster)
## No test: 
multiHeatMatrix(sml.ordered)
## End(No test)



cleanEx()
nameEx("plotGeneAnnotation-methods")
### * plotGeneAnnotation-methods

flush(stderr()); flush(stdout())

### Name: plotGeneAnnotation
### Title: Plots the enrichment of each feature in the set in the gene
###   annotation
### Aliases: plotGeneAnnotation
###   plotGeneAnnotation,AnnotationByFeature-method
###   plotGeneAnnotation,list-method

### ** Examples

library(GenomicRanges)
data(cage)
data(cpgi)

cage$tpm = NULL
gl = GRangesList(cage=cage, cpgi=cpgi)

bed.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
annot = annotateWithGeneParts(gl, gene.parts, intersect.chr=TRUE)
## No test: 
plotGeneAnnotation(annot)
## End(No test)



cleanEx()
nameEx("plotMeta")
### * plotMeta

flush(stderr()); flush(stdout())

### Name: plotMeta
### Title: Line plot(s) for meta-region profiles
### Aliases: plotMeta

### ** Examples

data(cage)
 data(promoters)
 scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE)

 data(cpgi)
 scores2=ScoreMatrix(target=cpgi,windows=promoters,strand.aware=TRUE)

 # create a new ScoreMatrixList
 x=new("ScoreMatrixList",list(scores1,scores2))
 ## No test: 
 plotMeta(mat=x,overlay=TRUE,main="my plotowski")

 # plot dispersion nd smooth central tendency and variation interval bands
 plotMeta(mat=x, centralTend="mean", dispersion="se", winsorize=c(0,99),
         main="Dispersion as interquartile band", lwd=4,
         smoothfun=function(x) stats::lowess(x, f = 1/5))
 
## End(No test)



cleanEx()
nameEx("plotTargetAnnotation-methods")
### * plotTargetAnnotation-methods

flush(stderr()); flush(stdout())

### Name: plotTargetAnnotation
### Title: Plot annotation categories from AnnotationByGeneParts or
###   AnnotationByFeature
### Aliases: plotTargetAnnotation
###   plotTargetAnnotation,AnnotationByFeature-method

### ** Examples

data(cage)

bed.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
gene.parts = readTranscriptFeatures(bed.file)
annot = annotateWithGeneParts(cage, gene.parts, intersect.chr=TRUE)
## No test: 
plotTargetAnnotation(annot)
## End(No test)



cleanEx()
nameEx("readBed")
### * readBed

flush(stderr()); flush(stdout())

### Name: readBed
### Title: Read a BED file and convert it to GRanges.
### Aliases: readBed

### ** Examples

my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
refseq = readBed(my.file,track.line=FALSE,remove.unusual=FALSE)
head(refseq)



cleanEx()
nameEx("readBroadPeak")
### * readBroadPeak

flush(stderr()); flush(stdout())

### Name: readBroadPeak
### Title: A function to read the Encode formatted broad peak file into a
###   GRanges object
### Aliases: readBroadPeak

### ** Examples

broad.peak.file = system.file('extdata',"ex.broadPeak", package='genomation')

broad.peak = readBroadPeak(broad.peak.file)
head(broad.peak)



cleanEx()
nameEx("readFeatureFlank-methods")
### * readFeatureFlank-methods

flush(stderr()); flush(stdout())

### Name: readFeatureFlank
### Title: A function to read-in genomic features and their upstream and
###   downstream adjecent regions such as CpG islands and their shores
### Aliases: readFeatureFlank readFeatureFlank,character-method

### ** Examples

cgi.path = system.file('extdata/chr21.CpGi.hg19.bed', package='genomation')
cgi.shores = readFeatureFlank(cgi.path)
cgi.shores



cleanEx()
nameEx("readGeneric")
### * readGeneric

flush(stderr()); flush(stdout())

### Name: readGeneric
### Title: Read a tabular file and convert it to GRanges.
### Aliases: readGeneric

### ** Examples

my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
refseq = readGeneric(my.file,chr=1,start=2,end=3,strand=NULL,
                      meta.cols=list(score=5,name=4),
                      keep.all.metadata=FALSE, zero.based=TRUE)
head(refseq)



cleanEx()
nameEx("readNarrowPeak")
### * readNarrowPeak

flush(stderr()); flush(stdout())

### Name: readNarrowPeak
### Title: A function to read the Encode formatted narrowPeak file into a
###   GRanges object
### Aliases: readNarrowPeak

### ** Examples

narrow.peak.file = system.file('extdata',"ex.narrowPeak", package='genomation')

narrow.peak = readBroadPeak(narrow.peak.file)
head(narrow.peak)



cleanEx()
nameEx("readTranscriptFeatures-methods")
### * readTranscriptFeatures-methods

flush(stderr()); flush(stdout())

### Name: readTranscriptFeatures
### Title: Function for reading exon intron and promoter structure from a
###   given bed file
### Aliases: readTranscriptFeatures readTranscriptFeatures,character-method

### ** Examples

my.bed12.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
my.bed12.file
feats = readTranscriptFeatures(my.bed12.file)
names(feats)
sapply(feats, head)



cleanEx()
nameEx("scaleScoreMatrix-methods")
### * scaleScoreMatrix-methods

flush(stderr()); flush(stdout())

### Name: scaleScoreMatrix
### Title: Scales the values in the matrix by rows and/or columns
### Aliases: scaleScoreMatrix scaleScoreMatrix,ScoreMatrix-method

### ** Examples

# scale the rows of a scoreMatrix object
library(GenomicRanges)
target = GRanges(rep(c(1,2),each=7), IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
         weight = rep(c(1,2),each=7),
         strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))
windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5),
           strand=c('-','+','-','+'))
sm = ScoreMatrix(target, windows)
ssm = scaleScoreMatrix(sm, rows=TRUE)



cleanEx()
nameEx("scaleScoreMatrixList")
### * scaleScoreMatrixList

flush(stderr()); flush(stdout())

### Name: scaleScoreMatrixList
### Title: Scale the ScoreMatrixList
### Aliases: scaleScoreMatrixList
###   scaleScoreMatrixList,ScoreMatrixList-method

### ** Examples

library(GenomicRanges)
data(cage)
data(cpgi)
data(promoters)

cage$tpm = NULL
targets = GRangesList(cage=cage, cpgi=cpgi)
sml = ScoreMatrixList(targets, promoters, bin.num=10, strand.aware=TRUE)
sml.scaled = scaleScoreMatrixList(sml, rows=TRUE)
sml.scaled
## No test: 
multiHeatMatrix(sml)
## End(No test)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
