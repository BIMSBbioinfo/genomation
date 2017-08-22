# test C++ binSummary function

# to do: 

library(Rcpp)
library(genomation)
library(rtracklayer)
sourceCpp("binSum.cpp")

# get data from compgenomr.github.com/book
bwFile="wgEncodeHaibTfbsA549.chr21.bw"

# read refseq file
ref= readTranscriptFeatures("refseq.hg19.chr21.bed") 
promoter.gr=ref$promoters




myFunc<-function(bwFile,promoter.gr){
  cov.bw=import(bwFile, which=promoter.gr,as = "RleList");
  myViews=Views(cov.bw,as(promoter.gr,"RangesList")); # get subsets of coverage
  
  mat = lapply(myViews,function(x) as.list((viewApply(x,as.vector,
                                                      simplify = FALSE))) )
  listSliceMean(do.call("c",mat),10);
}

bigp=rep(promoter.gr,20) # make a realistic number of windows 

library(microbenchmark)

microbenchmark(
  myFunc(bwFile,bigp),
  ScoreMatrixBin(bwFile, bigp,bin.num = 10,type="bigWig")
)

