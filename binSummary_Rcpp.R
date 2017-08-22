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
  mcols(promoter.gr)$X_rank = 1:length(promoter.gr);
  cov.bw=import(bwFile, which=promoter.gr,as = "RleList");
  myViews=Views(cov.bw,as(promoter.gr,"RangesList")); # get subsets of coverage
  
  Xrank <- myViews$chr21@elementMetadata$X_rank
  
  mat = lapply(myViews,function(x) as.list((viewApply(x,as.vector,
                                                      simplify = FALSE))) )
  
  mat_res = listSliceMean(do.call("c",mat),10);
  ## copied from scoreMatrix.R
  # get the ranks of windows, when things are reorganized by as(...,"RangesList")
  r.list=split(mcols(promoter.gr)[,"X_rank"], as.vector(seqnames(promoter.gr))  )
  r.list=r.list[order(names(r.list))]
  ranks=do.call("c",r.list)
  rownames(mat_res) = ranks
  mat_res = mat_res[order(ranks),] 
  return(mat_res)

}

myFunc2<-function(bwFile,promoter.gr){
  cov.bw=import(bwFile, which=promoter.gr,as = "RleList");
  mcols(promoter.gr)$X_rank = 1:length(promoter.gr)
  
  myViews=Views(cov.bw,as(promoter.gr,"RangesList")); # get subsets of coverage
  
  Xrank <- myViews$chr21@elementMetadata$X_rank
  mat = lapply(myViews,function(x) as.list((viewApply(x,as.vector,
                                                      simplify = FALSE))) )
  mat_res = listSliceMean2(do.call("c",mat),10);
  
  ## copied from scoreMatrix.R
  # get the ranks of windows, when things are reorganized by as(...,"RangesList")
  r.list=split(mcols(promoter.gr)[,"X_rank"], as.vector(seqnames(promoter.gr))  )
  r.list=r.list[order(names(r.list))]
  ranks=do.call("c",r.list)
  rownames(mat_res) = ranks
  mat_res = mat_res[order(ranks),] 
  return(mat_res)
}

bigp=rep(promoter.gr,20) # make a realistic number of windows 

library(microbenchmark)

microbenchmark(
  myFunc(bwFile,bigp),
  ScoreMatrixBin(bwFile, bigp,bin.num = 10,type="bigWig"),
  myFunc2(bwFile,bigp)
)

