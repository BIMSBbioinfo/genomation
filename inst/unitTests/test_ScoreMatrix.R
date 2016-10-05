# ---------------------------------------------------------------------------- #
# test for ScoreMatrix function
test_ScoreMatrix_RleList_GRanges = function()
	{

		# input RleList
		rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		
		#1. test for proper workings
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)), strand=c('+','-','+','-'))
		m1 = as(matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T), 'ScoreMatrix')
		rownames(m1) = 1:4
		checkIdentical(ScoreMatrix(rl, gr1, strand.aware = FALSE), m1)
		
		m2 = as(matrix(c(1,1,1,3,2,2,4,4,4,6,5,5), ncol=3, byrow=T), 'ScoreMatrix')
		rownames(m2) = 1:4
		checkIdentical(ScoreMatrix(rl, gr1, strand.aware = TRUE), m2)
		
		#2. test for different lengths
		gr2 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,9)))
		checkException(ScoreMatrix(rl, gr2), silent=TRUE)
    #"width of 'windows' are not equal, provide 'windows' with equal widths")
		
		#3. test for removing outliers
		gr3 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,8,1),c(3,7,10,3)))
		m3 = as(matrix(c(1,1,1,2,2,3,4,4,4), ncol=3, byrow=T), 'ScoreMatrix')
		rownames(m3) = c(1,2,4)
		checkIdentical(ScoreMatrix(rl, gr3), m3)
		
		#4. test for different chromosomes
		gr4 = GRanges(rep(c('chr1','chr3'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		m4 = as(matrix(c(1,1,1,2,2,3), ncol=3, byrow=T), 'ScoreMatrix')
		rownames(m4) = 1:2
		checkIdentical(ScoreMatrix(rl, gr4, strand.aware = FALSE), m4)
		
		#5. no overlapping chromosomes
		gr5 = GRanges(rep(c('chr3','chr4'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		checkException(ScoreMatrix(rl, gr5), silent=TRUE)
    # "All windows fell have coordinates outside chromosome boundaries"
		
		#6. only windows that fall off chromosomes
		gr6 = GRanges(rep(c('chr1','chr4'), each=2), IRanges(c(-1,8,1,5),c(1,10,3,7)))
		checkException(ScoreMatrix(rl, gr6), silent=TRUE)
    # "All windows fell have coordinates outside chromosome boundaries")
		
	}

# ---------------------------------------------------------------------------- #
test_ScoreMatrix_GRanges_GRanges = function()
{
  target = GRanges(rep(c(1,2),each=6), 
                   IRanges(rep(c(1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=6))
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
          
  # -----------------------------------------------#
  # usage
  # normal function works
  s1 = ScoreMatrix(target=target, windows=windows)
  m1 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  rownames(m1) = 1:4
  m1 = as(m1, 'ScoreMatrix')
  checkEquals(s1,m1)
          
  # function with weight col works
  s2 = ScoreMatrix(target=target, windows=windows, weight.col='weight')
  m2 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  m2[3:4,] = m2[3:4,]*2
  rownames(m2) = 1:4
  m2 = as(m2, 'ScoreMatrix')
  checkEquals(s2,m2)  
          
  #strand aware
  s3 = ScoreMatrix(target=target, windows=windows, strand.aware=T)
  m3 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  rownames(m3) = 1:4
  m3[c(1,3),] = rev(m3[c(1,3),])
  m3 = as(m3, 'ScoreMatrix')
  checkEquals(s3,m3)
          
  # -----------------------------------------------#
  # errors
  checkException(ScoreMatrix(target, windows, weight.col=''), silent=TRUE)
          
}


# ---------------------------------------------------------------------------- #
test_ScoreMatrix_character_GRanges = function()
{
  target = GRanges(rep(c(1,2),each=7), 
                   IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=7), 
                   strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))  
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
  
  target.paired.end = GRanges(rep(1,each=12), 
                              IRanges(start=c(1,1,2,2,3,3,7,7,8,8,9,9), end=c(19, 19, 20, 20, 21, 21,22,22,23,23,24,24)),
                              strand=rep("+", times=12))
  windows.paired.end = GRanges(rep(c(1),each=4), IRanges(c(7,8,20, 21), width=10), 
                               strand=c('+','+','+','+'))
    

  # -----------------------------------------------#
  # usage
  # bam file
  bam.file = system.file('unitTests/test.bam', package='genomation')
  s1 = ScoreMatrix(bam.file, windows, type='bam')
  m1 = ScoreMatrix(target, windows)
  checkEquals(s1,m1)
  
  # bam file, rpm=TRUE
  s2 = ScoreMatrix(bam.file, windows, type='bam', rpm=TRUE)
  tot = 1e6/sum(idxStats(normalizePath(bam.file))$mapped)
  m2 = m1*tot
  checkEquals(s2, m2)
  
  #bam file, rpm=FALSE, unique=TRUE
  s3 = ScoreMatrix(bam.file, windows, type='bam', unique=TRUE)
  m3 = ScoreMatrix(unique(target), windows)
  checkEquals(s3,m3)
  
  #bam file, rpm=FALSE, unique=TRUE, extend=1
  s4 = ScoreMatrix(bam.file, windows, type='bam', unique=TRUE, extend=1)
  m4 = ScoreMatrix(resize(unique(target), width=1), windows)
  checkEquals(s4,m4)
  
  # bam file with paired-end reads
  bam.pe.file = system.file('unitTests/test_pairedend.bam', package='genomation')
  s5 = ScoreMatrix(bam.pe.file, windows.paired.end, type='bam', bam.paired.end=TRUE)
  names(target.paired.end) <-c("1:1-19","1:1-19.1","1:2-20","1:2-20","1:3-21","1:3-21",
                               "1:7-22","1:7-22","1:8-23","1:8-23","1:9-24","1:9-24")
  target.paired.end = target.paired.end[!duplicated(names(target.paired.end))]
  m5 = ScoreMatrix(target.paired.end, windows.paired.end, bam.paired.end=TRUE)
  checkEquals(s5,m5)
  
  # test library.size argument
  s6 = ScoreMatrix(bam.pe.file, windows.paired.end, type='bam', bam.paired.end=TRUE,
                   rpm=TRUE, library.size=14)
  s7 = ScoreMatrix(bam.pe.file, windows.paired.end, type='bam', bam.paired.end=TRUE,
                   rpm=TRUE)
  checkEquals(s6,s7)

  # -----------------------------------------------#
  # errors
  # error upon not specifying the file
  checkException(ScoreMatrix('',windows), silent=TRUE)
  
  # error upon not specifying the format
  checkException(ScoreMatrix(bam.file, target), silent=TRUE)
  
}

# ---------------------------------------------------------------------------- #
test_ScoreMatrix_character_GRanges_bigWig = function()
{
  if (.Platform$OS.type == "windows")
    return()
  test_bw <- system.file("unitTests/test.bw", package = "genomation")
          
  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
  s = ScoreMatrix(test_bw, g, type='bigWig')
          
  m = matrix(-1, ncol=10, nrow=6)
  m[6,-1] = -0.75
  rownames(m) = 1:6
  m = as(m, 'ScoreMatrix')
  checkEquals(s, m)
}


# ---------------------------------------------------------------------------- #
# test arithmetic and logic operations on ScoreMatrix
test_ScoreMatrix_Ops = function()
{
  target1 = GRanges(rep(c(1,2),times=6), 
                   IRanges(rep(c(1,2,3,7,8,9), times=2), width=5),
                   weight = rep(c(1,2),each=6))
  target2 = GRanges(rep(c(1,2),times=12), 
                   IRanges(rep(c(1,2,3,7,8,9), times=4), width=5),
                   weight = rep(c(1,2),each=12))
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
                    strand=c('-','+','-','+'))
  
  s1 = ScoreMatrix(target=target1, windows=windows)
  s2 = ScoreMatrix(target=target2, windows=windows)
  checkEquals(s1*2,s2)

  checkEquals(sum(s1==s1), nrow(s1)*ncol(s1))
}
	
# ---------------------------------------------------------------------------- #
# test for constrainRanges
# test_that("constrainRanges works",
# 	{
# 		target = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
# 		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,8),c(3,7,3,10)))
# 		gt1 = GRanges(c('chr1','chr1','chr2'), IRanges(c(1,5,1),c(3,7,3)))
# 		expect_identical(constrainRanges(target, gr1), gt1)
# 	}
# )

# ---------------------------------------------------------------------------- #
# test for binMatrix
test_binMatrix = function()
{
	m1 = as(matrix(rep(1:4, 4), ncol=4), 'ScoreMatrix')
	
	# no binning
	checkEquals(binMatrix(m1), m1)
	
	# nbins 2, default function
	mb1 = as(matrix(c(1,1,2,2,3,3,4,4), ncol=2, byrow=T), 'ScoreMatrix')
	checkEquals(binMatrix(m1, bin.num=2), mb1)
	
	# nbins 5
	checkException(binMatrix(m1, bin.num=5), silent=TRUE)
	
	# nbins 2, not default function
	m2 = as(matrix(rep(1:4, 4), ncol=4, byrow=T), 'ScoreMatrix')
	mb2 = as(matrix(rep(c(2,4), 4), ncol=2, byrow=T), 'ScoreMatrix')
	checkEquals(binMatrix(m2, bin.num=2, max), mb2)	
}


# ---------------------------------------------------------------------------- #
test_scaleScoreMatrix = function()
{
  s = as(matrix(1:9, ncol=3),'ScoreMatrix')
          
  s1=scaleScoreMatrix(s, rows=T, columns=F)
  m1 = as(matrix(c(rep(-0.4285714,3),rep(0,3), rep(0.4285714,3)), 
                 nrow=3),'ScoreMatrix')
  checkEquals(s1,m1, tolerance=1e-5)
          
  s2=scaleScoreMatrix(s, rows=F, columns=T)
  m2 = as(matrix(c(rep(-0.3333333,3),rep(0,3), rep(0.3333333,3)), 
                 nrow=3, byrow=T),'ScoreMatrix')
  checkEquals(s2,m2, tolerance=1e-5)
          
  s3=scaleScoreMatrix(s, rows=T, columns=T)
  m3 = as(matrix(0,nrow=3,ncol=3),'ScoreMatrix')
  checkEquals(s3,m3, tolerance=1e-5)
}