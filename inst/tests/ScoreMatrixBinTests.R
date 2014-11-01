# ---------------------------------------------------------------------------- #
# test for ScoreMatrix function
test_that("ScoreMatrixBin: RleList, GRanges works",
{
  
  # -----------------------------------------------#
  # usage
  # input RleList
  rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
  
  #1. test for proper workings
  gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)), 
                strand=c('+','-','+','-'))
  m1 = matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T)
  m1 = as(cbind(rowMeans(m1[,1:2]), rowMeans(m1[,2:3])),'ScoreMatrix')
  rownames(m1) = 1:4
  s1 = ScoreMatrixBin(rl, gr1, bin.num=2)
  expect_equal(s1, m1)
  
  #2. test for different bin.op
  for(fun in c('min', 'max', 'median')){
    message('testing:', fun)
    m2 = matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T)
    m2 = as(cbind(apply(m2[,1:2],1,match.fun(fun)), 
                  apply(m2[,2:3],1,match.fun(fun))),'ScoreMatrix')
    rownames(m2) = 1:4
    s2 = ScoreMatrixBin(rl, gr1, bin.num=2, bin.op=fun)
    expect_equal(s2, m2)  
  }
  
  #3. test strand aware
  m3 = matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T)
  m3 = as(cbind(rowMeans(m3[,1:2]), rowMeans(m3[,2:3])),'ScoreMatrix')
  rownames(m3) = 1:4
  m3[c(2,4),] = m3[c(2,4),2:1]
  s3 = ScoreMatrixBin(rl, gr1, bin.num=2, strand.aware=T)
  expect_equal(s3, m3)
  
  # -----------------------------------------------#
  # errors
  # error for removing all bins
  expect_error(SCoreMatrixBin(rl, gr1, bin.num=2))
  
  gr5 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(5,6,5,6)), 
                strand=c('+','-','+','-'))
  expect_warning(ScoreMatrixBin(rl, gr5, bin.num=3))
  
})

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
  bam.file = system.file('tests/test.bam', package='genomation')
  s1 = ScoreMatrixBin(bam.file, windows, type='bam', bin.num=2)
  m1 = ScoreMatrixBin(target, windows, bin.num=2)
  expect_equal(s1,m1)
  
  # bam file, rpm=TRUE
  s2 = ScoreMatrixBin(bam.file, windows,bin.num=2, type='bam', rpm=TRUE)
  tot = 1e6/countBam(BamFile(bam.file))$records
  m2 = m1*tot
  expect_equal(s2, m2)
  
  #bam file, rpm=FALSE, unique=TRUE
  s3 = ScoreMatrixBin(bam.file, windows,bin.num=2, type='bam', unique=T)
  m3 = ScoreMatrixBin(unique(target), windows, bin.num=2)
  expect_equal(s3,m3)
  
  #bam file, rpm=FALSE, unique=TRUE, extend=1
  s4 = ScoreMatrixBin(bam.file, windows, type='bam',bin.num=2, unique=T, extend=1)
  m4 = ScoreMatrixBin(resize(unique(target), width=1), windows, bin.num=2)
  expect_equal(s4,m4)
  
  #bigWig file - missing
  
  # -----------------------------------------------#
  # errors
  # error upon not specifying the file
  expect_error(ScoreMatrixBin('',windows))
  
  # error upon not specifying the format
  expect_error(ScoreMatrixBin(bam.file, target))
})

# ---------------------------------------------------------------------------- #
test_that("ScoreMatrix:character, GRanges, type='bigWig' works".
{
  test_path <- system.file("tests/test.bw", package = "genomation")
  b = import(test_bw, asRangedData=F)
  
  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=5))
  s = ScoreMatrixBin(test_bw, g, type='bigWig', bin.num=2)
  
  m = matrix(-1, ncol=2, nrow=6)
  m[6,1] = -0.833333333333
  m[6,2] = -0.75
  rownames(m) = 1:6
  m = as(m, 'ScoreMatrix')
  expect_equal(s, m)
})