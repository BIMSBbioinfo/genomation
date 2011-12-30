#######################
# tests and development scripts to be deleted before release
# ONLY for developer use!!!!!!!!
#######################

# get data
#save(slentghs, feature, file = "~/Dropbox/PAPERS/R-devel/genomation_trials/feature.RData")
#load("~/Dropbox/PAPERS/R-devel/genomation_trials/feature.RData")
#feat=GRanges(seqnames=as.character( seqnames(feature)),
#                                        ranges  =ranges( feature),
#                                        strand  =strand(feature),
#                                        elementMetadata(feature)
#            )
#feature=feat
#chr.len=slentghs

#library(methylKit)
#chr.gaps=read.bed("~/Downloads/gaps.bed")

#cage=read.bed("~/Dropbox/genomation_data/L3_clusters.humanFantom4.tagTh2.cs.tagcnt.bed")

#detach(package:methylKit)
#save(cage,chr.gaps,chr.len, feature, file = "~/Dropbox/PAPERS/R-devel/genomation_trials/feature.RData")

load("~/Dropbox/PAPERS/R-devel/genomation_trials/data/feature.RData")

library(devtools)
load_all("~/Dropbox/PAPERS/R-devel/genomation_trials")
#source("~/Dropbox/PAPERS/R-devel/genomation_trials/R/getRandomEnrichment.R")
#tester and development script

randomize.feature(feature)
randomize.feature(feature,chrom.sizes=chr.len)
randomize.feature(feature,chrom.sizes=chr.len,exclude=chr.gaps)


trace("randomize.feature", browser, exit=browser, signature = c("GRanges"))


trace("getRandomEnrichment", browser, exit=browser, signature = c("GRanges","GRanges"))
getRandomEnrichment(target=cage,query=feature,randomizations=10)
system.time(getRandomEnrichment(target=cage,query=feature,randomizations=10)) # 11.2 seconds for my macbook

x=modCoverage(cage,col.name="score",multiply=1000,add=1)
is(x,"RleList") # TRUE
is(x,"modRleList") # TRUE 
