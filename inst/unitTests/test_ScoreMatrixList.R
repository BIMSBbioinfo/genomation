# ---------------------------------------------------------------------------- #
# scoreMatrixList: standard
test_scoreMatrixList = function()
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
  checkEquals(sml1, tml1)
  
  #2. ScoreMatrixList: list, GRanges
  sml2 = ScoreMatrixList(tar.list, target)
  tml2 = as(list(tar1=ScoreMatrix(tar.list$tar1, target), 
                 tar2=ScoreMatrix(tar.list$tar1, target)),
            'ScoreMatrixList')
  checkEquals(sml2, tml2)
  
  #3. ScoreMatrixList: GRangesList, GRanges
  sml3 = ScoreMatrixList(tar.glist, target)
  tml3 = as(list(tar1=ScoreMatrix(tar.list$tar1, target), 
                 tar2=ScoreMatrix(tar.list$tar1, target)),
            'ScoreMatrixList')
  checkEquals(sml3, tml3)
  
  #4. ScoreMatrixList: list, GRanges, weight.col
  sml4 = ScoreMatrixList(tar.list, target, weight.col='weight')
  tml4 = as(list(tar1=ScoreMatrix(tar.list$tar1, target, weight.col='weight'), 
                 tar2=ScoreMatrix(tar.list$tar1, target, weight.col='weight')),
            'ScoreMatrixList')
  checkEquals(sml4, tml4)
  
  #5. ScoreMatrixList: list, GRanges, strand.aware
  sml5 = ScoreMatrixList(tar.list, target, strand.aware=TRUE)
  tml5 = as(list(tar1=ScoreMatrix(tar.list$tar1, target, strand.aware=TRUE), 
                 tar2=ScoreMatrix(tar.list$tar1, target, strand.aware=TRUE)),
            'ScoreMatrixList')
  checkEquals(sml5, tml5)
  
  #6. ScoreMatrixList: list, GRanges, Bin
  sml6 = ScoreMatrixList(tar.list, target, bin.num=2)
  tml6 = as(list(tar1=ScoreMatrixBin(tar.list$tar1, target, bin.num=2), 
                 tar2=ScoreMatrixBin(tar.list$tar1, target, bin.num=2)),
            'ScoreMatrixList')
  checkEquals(sml6, tml6)
  
  #7. ScoreMatrixList: list, GRanges, Bin, bin.op='min'
  sml7 = ScoreMatrixList(tar.list, target, bin.num=2, bin.op='min')
  tml7 = as(list(tar1=ScoreMatrixBin(tar.list$tar1, target, bin.num=2, bin.op='min'), 
                 tar2=ScoreMatrixBin(tar.list$tar1, target, bin.num=2, bin.op='min')),
            'ScoreMatrixList')
  checkEquals(sml7, tml7)
  
  
  
  # -----------------------------------------------#
  #1. error given a single granges object
  checkException(ScoreMatrixList(target, windows), silent=TRUE)          
  
  #2. error for empty target
  checkException(ScoreMatrixList(NULL, windows), silent=TRUE)
  
  #3. error for empty windows
  checkException(ScoreMatrixList(tar.list, NULL), silent=TRUE)
  
  #4. error for character not a file
  checkException(ScoreMatrixList(c('',''), windows), silent=TRUE) 
}

# ---------------------------------------------------------------------------- #
# scoreMatrixList: character, bigWig
test_ScoreMatrixList_BigWig = function()
{
  if (.Platform$OS.type == "windows")
    return()
  test_bw <- system.file("unitTests/test.bw", package = "genomation")
          
  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
  s = ScoreMatrixList(c(test_bw, test_bw), g, type='bigWig')
          
  m = matrix(-1, ncol=10, nrow=6)
  m[6,-1] = -0.75
  rownames(m) = 1:6
  m = as(m, 'ScoreMatrix')
  l = as(list(test.bw=m, test.bw=m), 'ScoreMatrixList')
  checkEquals(s, l)
}

# ---------------------------------------------------------------------------- #
# scoreMatrixList: character, bigWig, bin
test_ScoreMatrixListBin_bigWig = function()
{
  if (.Platform$OS.type == "windows")
    return()
  test_bw <- system.file("unitTests/test.bw", package = "genomation")
          
  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
#   b = import(test_bw, which=g)
#   covs = coverage(b, weight=b$score)        
#   s = ScoreMatrixBin(covs, g, bin.num=5)
          
  s = ScoreMatrixList(c(test_bw, test_bw), g,bin.num=5, type='bigWig')
          
  m = matrix(-1, ncol=10, nrow=6)
  m[6,-1] = -0.75
  rownames(m) = 1:6
  m = as(m, 'ScoreMatrix')
  m = binMatrix(m, bin.num=5)
  l = as(list(test.bw=m, test.bw=m), 'ScoreMatrixList')
  checkEquals(s, l)
}

# ---------------------------------------------------------------------------- #
test_ScoreMatrixList_Bam = function()
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
  bam.file = system.file('unitTests/test.bam', package='genomation')
  bam.files = c(bam.file, bam.file)
  sml1 = ScoreMatrixList(bam.files, windows, type='bam')
  tml1 = ScoreMatrixList(tar.list, windows)
  names(tml1) = basename(bam.files)
  checkEquals(sml1,tml1)
          
  # bam file, rpm=TRUE
  s2 = ScoreMatrixList(bam.files, windows, type='bam', rpm=TRUE)
  tot= 1e6/countBam(BamFile(bam.file))$records
  m2 = ScoreMatrixList(lapply(sml1, function(x)x*tot))
  checkEquals(s2, m2)
          
  #bam file, rpm=FALSE, unique=TRUE
  s3 = ScoreMatrix(bam.file, windows, type='bam', unique=T)
  m3 = ScoreMatrix(unique(target), windows)
  checkEquals(s3,m3)
          
  #bam file, rpm=FALSE, unique=TRUE, extend=1
  s4 = ScoreMatrix(bam.file, windows, type='bam', unique=T, extend=1)
  m4 = ScoreMatrix(resize(unique(target), width=1), windows)
  checkEquals(s4,m4)    
  
  # bam file, rpm=TRUE, library.size
  s3 = ScoreMatrixList(bam.files, windows, type='bam', rpm=TRUE)
  tot= sum(countBam(BamFile(bam.file))$records)
  s4 = ScoreMatrixList(bam.files, windows, type='bam', rpm=TRUE,
                       library.size=c(tot,tot))
  checkEquals(s3, s4)
  
}

# ---------------------------------------------------------------------------- #
# test arithmetic on ScoreMatrixList
test_ScoreMatrixList_Ops = function()
{
  target = GRanges(rep(c(1,2),each=7), 
                   IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=7), 
                   strand=c('-', '-', '-', '-', '+', '-', '+', '-', 
                            '-', '-', '-', '-', '-', '+'))
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
  
  tar.list = list(tar1=target, tar2=target)
  sml = ScoreMatrixList(tar.list, windows)
  sml1 = sml*10
  sml2 = as(lapply(sml, function(x) x*10), "ScoreMatrixList")
  checkEquals(sml1, sml2, checkNames=FALSE)
}




