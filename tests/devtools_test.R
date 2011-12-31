#######################
# tests and development scripts to be deleted before release
# ONLY for developer use!!!!!!!!
#######################
library(roxygen2)
roxygenise("~/Dropbox/PAPERS/R-devel/genomation",roclets = c("rd"))

load("~/Dropbox/PAPERS/R-devel/genomation/data/feature.RData")
# cage: cage TCs hg18
# feature: CpG islands hg18
# chr.len: chromosome lengths hg18
# chr.gaps: assembly gaps on hg18



library(devtools)
load_all("~/Dropbox/PAPERS/R-devel/genomation")




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
