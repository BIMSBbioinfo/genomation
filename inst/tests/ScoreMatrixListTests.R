# ---------------------------------------------------------------------------- #
# scoreMatrixList: standard
test_that("scoreMatrixList basic works"),
{
  target = GRanges(rep(c(1,2),each=7), 
                   IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=7), 
                   strand=c('-', '-', '-', '-', '+', '-', '+', '-', 
                            '-', '-', '-', '-', '-', '+'))
  tar.list = list(tar1=target, tar2=target)
  tar.glist = GRangesList(tar1=target, tar2=target)
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
  
  # -----------------------------------------------#
  # usage
  #1. return sml
  m = as(matrix(0, ncol=5, nrow=4), 'ScoreMatrix')
  ml1 = list(m1=m, m2=m)
  sml1 = ScoreMatrixList(ml1)
  tml1 = as(ml1, 'ScoreMatrixList')
  expect_equal(sml1, tml1)
  
  #2. ScoreMatrixList: list, GRanges
  sml2 = ScoreMatrixList(tar.list, target)
  tml2 = as(list(tar1=ScoreMatrix(tar.list$tar1, target), 
                 tar2=ScoreMatrix(tar.list$tar1, target)),
            'ScoreMatrixList')
  expect_equal(sml2, tml2)
  
  #3. ScoreMatrixList: GRangesList, GRanges
  sml3 = ScoreMatrixList(tar.glist, target)
  tml3 = as(list(tar1=ScoreMatrix(tar.list$tar1, target), 
                 tar2=ScoreMatrix(tar.list$tar1, target)),
            'ScoreMatrixList')
  expect_equal(sml3, tml3)
  
  #4. ScoreMatrixList: list, GRanges, weight.col
  sml4 = ScoreMatrixList(tar.list, target, weight.col='weight')
  tml4 = as(list(tar1=ScoreMatrix(tar.list$tar1, target, weight.col='weight'), 
                 tar2=ScoreMatrix(tar.list$tar1, target, weight.col='weight')),
            'ScoreMatrixList')
  expect_equal(sml4, tml4)
  
  #5. ScoreMatrixList: list, GRanges, strand.aware
  sml5 = ScoreMatrixList(tar.list, target, strand.aware=TRUE)
  tml5 = as(list(tar1=ScoreMatrix(tar.list$tar1, target, strand.aware=TRUE), 
                 tar2=ScoreMatrix(tar.list$tar1, target, strand.aware=TRUE)),
            'ScoreMatrixList')
  expect_equal(sml5, tml5)
  
  #6. ScoreMatrixList: list, GRanges, Bin
  sml6 = ScoreMatrixList(tar.list, target, bin.num=2)
  tml6 = as(list(tar1=ScoreMatrixBin(tar.list$tar1, target, bin.num=2), 
                 tar2=ScoreMatrixBin(tar.list$tar1, target, bin.num=2)),
            'ScoreMatrixList')
  expect_equal(sml6, tml6)
  
  #7. ScoreMatrixList: list, GRanges, Bin, bin.op='min'
  sml7 = ScoreMatrixList(tar.list, target, bin.num=2, bin.op='min')
  tml7 = as(list(tar1=ScoreMatrixBin(tar.list$tar1, target, bin.num=2, bin.op='min'), 
                 tar2=ScoreMatrixBin(tar.list$tar1, target, bin.num=2, bin.op='min')),
            'ScoreMatrixList')
  expect_equal(sml7, tml7)
  
  
  
  # -----------------------------------------------#
  #1. error given a single granges object
  expect_error(ScoreMatrixList(target, windows))          
  
  #2. error for empty target
  expect_error(ScoreMatrixList(NULL, windows))
  
  #3. error for empty windows
  expect_error(ScoreMatrixList(tar.list, NULL))
  
  #4. error for character not a file
  expect_error(ScoreMatrixList(c('',''), windows))
  
})

# ---------------------------------------------------------------------------- #
# scoreMatrixList: character, bigWig
test_that("ScoreMatrixList BigWig works",
{
  test_bw <- system.file("tests/test.bw", package = "genomation")
          
  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
  s = ScoreMatrixList(c(test_bw, test_bw), g, type='bigWig')
          
  m = matrix(-1, ncol=10, nrow=6)
  m[6,-1] = -0.75
  rownames(m) = 1:6
  m = as(m, 'ScoreMatrix')
  l = as(list(test.bw=m, test.bw=m), 'ScoreMatrixList')
  expect_equal(s, l)
})

# ---------------------------------------------------------------------------- #
# scoreMatrixList: character, bigWig, bin
test_that("ScoreMatrixList Bin bigWig works",
{
  test_bw <- system.file("tests/test.bw", package = "genomation")
          
  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
  b = import(test_bw, asRangedData=FALSE, which=g)
  covs = coverage(b, weight=b$score)        
  s = ScoreMatrixBin(covs, g, bin.num=5)
          
  s = ScoreMatrixList(c(test_bw, test_bw), g,bin.num=5, type='bigWig')
          
  m = matrix(-1, ncol=10, nrow=6)
  m[6,-1] = -0.75
  rownames(m) = 1:6
  m = as(m, 'ScoreMatrix')
  m = binMatrix(m, bin.num=5)
  l = as(list(test.bw=m, test.bw=m), 'ScoreMatrixList')
  expect_equal(s, l)
})

# ---------------------------------------------------------------------------- #
test_that("ScoreMatrixList Bam works",
{
  target = GRanges(rep(c(1,2),each=7), 
                   IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
                           weight = rep(c(1,2),each=7), 
                           strand=c('-', '-', '-', '-', '+', '-', '+', '-', 
                                    '-', '-', '-', '-', '-', '+'))
  tar.list = list(tar1=target, tar2=target)
  tar.glist = GRangesList(tar1=target, tar2=target)
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
          
  # -----------------------------------------------#
  # usage
  # bam file
  bam.file = system.file('tests/test.bam', package='genomation')
  bam.files = c(bam.file, bam.file)
  sml1 = ScoreMatrixList(bam.files, windows, type='bam')
  tml1 = ScoreMatrixList(tar.list, windows)
  names(tml1) = basename(bam.files)
  expect_equal(sml1,tml1)
          
  # bam file, rpm=TRUE
  s2 = ScoreMatrixList(bam.files, windows, type='bam', rpm=TRUE)
  tot= 1e6/countBam(BamFile(bam.file))$records
  m2 = ScoreMatrixList(lapply(sml1, function(x)x*tot))
  expect_equal(s2, m2)
          
  #bam file, rpm=FALSE, unique=TRUE
  s3 = ScoreMatrix(bam.file, windows, type='bam', unique=T)
  m3 = ScoreMatrix(unique(target), windows)
  expect_equal(s3,m3)
          
  #bam file, rpm=FALSE, unique=TRUE, extend=1
  s4 = ScoreMatrix(bam.file, windows, type='bam', unique=T, extend=1)
  m4 = ScoreMatrix(resize(unique(target), width=1), windows)
  expect_equal(s4,m4)        
})

