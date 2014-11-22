
# 
library(GenomicRanges)


# GET CAGE DATA
cage=read.table("~/Dropbox/annotation_data/wgEncodeRikenCageNhekNucleusPapTssHmmV3.bedRnaElements.txt",
          header=FALSE,sep="\t",stringsAsFactors=FALSE)
cage=cage[cage[,1] %in% c("chr22","chr21"),]
cage=GRanges(seqnames=cage[,1],ranges=IRanges(cage[,2],cage[,3]),
                  strand=cage[,6],tpm=cage[,7] )

# GET REFSEQ GENES and produce GRanges
ref=read.table("~/Dropbox/annotation_data/hg19.refseq.bed",
           header=FALSE,sep="\t",stringsAsFactors=FALSE)
ref=ref[ref[,1] %in% c("chr22","chr21"),]

# get promoters
prom=ref
prom[prom[,6]=="+",3]=prom[prom[,6]=="+",2]
prom[prom[,6]=="-",2]=prom[prom[,6]=="-",3]
prom[,2]=prom[,2]-1000
prom[,3]=prom[,3]+1000
promoters=GRanges(seqnames=prom[,1],ranges=IRanges(prom[,2],prom[,3]),
                  strand=prom[,6] )
promoters=unique(promoters)

# get genes
genes=GRanges(seqnames=ref[,1],ranges=IRanges(ref[,2],ref[,3]),
                  strand=ref[,6],name=ref[,4])

# GET CPG ISLANDS and produce GRanges
cpg=read.table("~/Dropbox/annotation_data/hg19.CpGi.bed",
               header=FALSE,sep="\t",stringsAsFactors=FALSE)
cpg=cpg[cpg[,1] %in% c("chr22","chr21"),]
cpgi=GRanges(seqnames=cpg[,1],ranges=IRanges(cpg[,2],cpg[,3]) )

# GET CHROMOSOME SIZES
chr.sizes=read.table("~/Dropbox/annotation_data/hg19.chrom.sizes",
                    header=FALSE,sep="\t",stringsAsFactors=FALSE)


save(cage     ,file="data/cage.RData")
save(cpgi     ,file="data/cpgi.RData")
save(genes    ,file="data/genes.RData")
save(promoters,file="data/promoters.RData")


