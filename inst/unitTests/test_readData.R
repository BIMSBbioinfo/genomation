# test for readGeneric function
test_readGeneric = function()
{
  tab.test1=system.file('unitTests/tab.test1', package='genomation')
  
  #1. test regular bed input
  r1 = readGeneric(tab.test1)
  g1 = GRanges(c('chr1','chr1'), IRanges(c(1,5), c(10,15)))
  checkIdentical(g1, r1)
  
  #2. test strand
  r2 = readGeneric(tab.test1, strand=4)
  g2 = GRanges(c('chr1','chr1'), IRanges(c(1,5), c(10,15)), strand=c('+','-'))
  checkIdentical(g2, r2)
  
  #3. all metadata colums
  r3 = readGeneric(tab.test1, strand=4, keep.all.metadata=TRUE)
  g3 = GRanges(c('chr1','chr1'), 
               IRanges(c(1,5), c(10,15)), 
               strand=c('+','-'),
               V5=as.integer(c(15,20)),
               V6=as.integer(c(20,25)))
  checkIdentical(g3,r3)
  
  #4. all named metadata colums
  r4 = readGeneric(tab.test1, strand=4, meta.col=list(score1=5, score2=6))
  g4 =  GRanges(c('chr1','chr1'), 
                IRanges(c(1,5), c(10,15)), 
                strand=c('+','-'),
                score1=as.integer(c(15,20)),
                score2=as.integer(c(20,25)))
  checkIdentical(g4,r4)
  
  #5. selected metadata columns
  r5 = readGeneric(tab.test1, strand=4, meta.cols=list(score1=6))
  g5 =  GRanges(c('chr1','chr1'), 
                IRanges(c(1,5), c(10,15)), 
                strand=c('+','-'),
                score1=as.integer(c(20,25)))
  checkIdentical(g5,r5)
  
  #6. test whether it can read a file containing a header
  tab.test2=system.file('unitTests/tab.test2', package='genomation')
  r6 = readGeneric(tab.test2, header=TRUE)
  g6 = GRanges(c('chr1','chr1'), IRanges(c(1,5), c(10,15)))
  checkIdentical(g6,r6)
  
  #7. header=TRUE, keep.all.metadata=TRUE
  r7 = readGeneric(tab.test2,strand=4, header=TRUE, keep.all.metadata=TRUE)
  
  
  #8. test whether it can read a file with permutted columns
  tab.test3=system.file('unitTests/tab.test3', package='genomation')
  r7 = readGeneric(tab.test3, chr=5, start=3, end=4, strand=6, 
                   meta.col=c(score1=1, score2=2), header=TRUE)
  g7 =  GRanges(c('chr1','chr1'), 
                IRanges(c(1,5), c(10,15)), 
                strand=c('+','-'),
                score1=as.integer(c(15,20)),
                score2=as.integer(c(20,25)))
  checkIdentical(g7,r7)
  
  
  #9. test whether it can read a compressed file (.gz, .zip)
  tab.test3=system.file('unitTests/tab.test3', package='genomation')
  tab.test3.zip=system.file('unitTests/tab.test3.zip', package='genomation')

  tab.test3.gz <- paste(tab.test3, ".gz", sep="")
  gzf <- gzfile(tab.test3.gz, "w")
  write.table(read.table(tab.test3, header=TRUE), 
              file=gzf, 
              sep = "\t", row.names = FALSE)
  close(gzf)
  
  r7.gz = readGeneric(tab.test3.gz, chr=5, start=3, end=4, strand=6, 
                      meta.col=c(score1=1, score2=2), header=TRUE)  
  r7.zip = readGeneric(tab.test3.zip, chr=5, start=3, end=4, strand=6, 
                       meta.col=c(score1=1, score2=2), header=TRUE)        
                       
  g7 =  GRanges(c('chr1','chr1'), 
                IRanges(c(1,5), c(10,15)), 
                strand=c('+','-'),
                score1=as.integer(c(15,20)),
                score2=as.integer(c(20,25)))
  
  checkIdentical(g7, r7.gz)
  checkIdentical(g7, r7.zip)
  
  if (file.exists(tab.test3.gz)) file.remove(tab.test3.gz)
  
  #10. test a file with UCSC header
  tab.test4=system.file('unitTests/tab.test4', package='genomation')
  r8 = readGeneric(tab.test4, chr=1, start=2, end=3, strand=6, 
                   meta.col=c(score=5,
                              name=4, 
                              thickStart=7,  
                              thickEnd=8, 
                              itemRgb=9), 
                   header=FALSE, skip=3,
                   zero.based=TRUE)
  r9 = readBed(tab.test4, track.line=3)
  r10 = readBed(tab.test4, track.line="auto")
  
  g8 =  GRanges(c('chr7','chr7','chr7'), 
                IRanges(c(1,11,16), c(10,15,20)), 
                strand=c('+','+','+'),
                score=as.integer(c(0,0,0)),
                name=c("Pos1","Pos2","Pos3"),
                thickStart=as.integer(c(0,10,15)),
                thickEnd=as.integer(c(10,15,20)),
                itemRgb=c("255,0,0","255,0,0","255,0,0"))
  
  checkIdentical(g8, r8)
  checkIdentical(g8, r9)
  checkIdentical(g8, r10)
}

test_gffToGRanges = function()
{
  library(GenomicRanges)
  tab.test = system.file('unitTests/test.gtf', package='genomation')
  gff1 = gffToGRanges(tab.test)
  checkIdentical(length(gff1), 3L)
 
  gff2 = gffToGRanges(tab.test, filter='exon')
  checkIdentical(length(gff2), 1L)
  
  gff3 = gffToGRanges(tab.test, ensembl=TRUE)
  checkIdentical(as.character(seqlevels(gff3)), 'chr1')
}

