# ---------------------------------------------------------------------------- #
test_that("ScoreMatrixBin:GRanges, GRanges works",
{
  target = GRanges(rep(c(1,2),each=6), 
                   IRanges(rep(c(1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=6))
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
  
  # -----------------------------------------------#
  # usage
  # normal function works
  s1 = ScoreMatrixBin(target=target, windows=windows, bin.num=2)
  m1 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  m1 = as(cbind(rowMeans(m1[,1:3]), rowMeans(m1[,3:5])), 'ScoreMatrix')
  rownames(m1) = 1:4
  expect_equal(s1,m1)
  
  # function with weight col works
  s2 = ScoreMatrixBin(target=target, windows=windows,bin.num=2, weight.col='weight')
  m2 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  m2[3:4,] = m2[3:4,]*2
  m2 = as(cbind(rowMeans(m2[,1:3]), rowMeans(m2[,3:5])), 'ScoreMatrix')
  rownames(m2) = 1:4
  expect_equal(s2,m2)  
  
  #strand aware
  s3 = ScoreMatrixBin(target=target, windows=windows, bin.num=2, strand.aware=T)
  m3 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  rownames(m3) = 1:4
  m3[c(1,3),] = rev(m3[c(1,3),])
  m3 = as(cbind(rowMeans(m3[,1:3]), rowMeans(m3[,3:5])), 'ScoreMatrix')
  rownames(m3) = 1:4
  expect_equal(s3,m3)
  
  # -----------------------------------------------#
  # errors
  expect_error(ScoreMatrixBin(target, windows, weight.col=''))
  
  # number of bins > ncol 
  expect_warning(ScoreMatrixBin(target, windows, strand.aware = FALSE, bin.num=10))
              
  
})


# ---------------------------------------------------------------------------- #
test_that("ScoreMatrixBin:character, GRanges works",
{
  target = GRanges(rep(c(1,2),each=7), 
                   IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=7), 
                   strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
  
  # -----------------------------------------------#
  # usage
  # bam file
  bam.file = system.file('extdata/test.bam', package='genomation')
  s1 = ScoreMatrixBin(bam.file, windows, type='bam', bin.num=2)
  m1 = ScoreMatrixBin(target, windows, bin.num=2)
  expect_equal(s1,m1)
  
  # bam file, rpm=TRUE
  s2 = ScoreMatrixBin(bam.file, windows,bin.num=2, type='bam', rpm=TRUE)
  tot = 1e6/countBam(BamFile(bam.file))$records
  m2 = m1*tot
  expect_equal(s2, m2)
  
  #bam file, rpm=FALSE, unique=TRUE
  s3 = ScoreMatrixBin(bam.file, windows, type='bam', unique=T)
  m3 = ScoreMatrixBin(unique(target), windows)
  expect_equal(s3,m3)
  
  #bam file, rpm=FALSE, unique=TRUE, extend=1
  s4 = ScoreMatrixBin(bam.file, windows, type='bam', unique=T, extend=1)
  m4 = ScoreMatrixBin(resize(unique(target), width=1), windows)
  expect_equal(s4,m4)
  
  #bigWig file - does not work on windows
  #bw.file = system.file('extdata/test.bw', package='genomation')
  #s5 = ScoreMatrix(bw.file, windows, type='bigWig')
  #m5 = ScoreMatrix(target, windows)
  #expect_equal(s1,m1)
  # -----------------------------------------------#
  # errors
  # error upon not specifying the file
  expect_error(ScoreMatrixBin('',windows))
  
  # error upon not specifying the format
  expect_error(ScoreMatrixBin(bam.file, target))
  
  # error for removing all bins
  expect_error(SCoreMatrixBin(bam.file, windows, type='bam', rpm=TRUE))
  
})