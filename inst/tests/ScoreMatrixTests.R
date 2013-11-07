# ---------------------------------------------------------------------------- #
# test for scoreMatrix function
test_that("scoreMatrix works",
	{
	
		# input RleList
		rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		
		#1. test for proper workings
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)), strand=c('+','-','+','-'))
		m1 = as(matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m1) = 1:4
		expect_equal(scoreMatrix(rl, gr1, strand.aware = FALSE), m1)
		
		m2 = as(matrix(c(1,1,1,3,2,2,4,4,4,6,5,5), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m2) = 1:4
		expect_equal(scoreMatrix(rl, gr1, strand.aware = TRUE), m2)
		
		#2. test for different lengths
		gr2 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,9)))
		expect_error(scoreMatrix(rl, gr2), "width of 'windows' are not equal, provide 'windows' with equal widths")
		
		#3. test for removing outliers
		gr3 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,8,1),c(3,7,10,3)))
		m3 = as(matrix(c(1,1,1,2,2,3,4,4,4), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m3) = c(1,2,4)
		expect_equal(scoreMatrix(rl, gr3), m3)
		
		#4. test for different chromosomes
		gr4 = GRanges(rep(c('chr1','chr3'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		m4 = as(matrix(c(1,1,1,2,2,3), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m4) = 1:2
		expect_equal(scoreMatrix(rl, gr4, strand.aware = FALSE), m4)
		
		#5. no overlapping chromosomes
		gr5 = GRanges(rep(c('chr3','chr4'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		expect_error(scoreMatrix(rl, gr5), "All windows fell have coordinates outside chromosome boundaries")
		
		#6. only windows that fall off chromosomes
		gr6 = GRanges(rep(c('chr1','chr4'), each=2), IRanges(c(-1,8,1,5),c(1,10,3,7)))
		expect_error(scoreMatrix(rl, gr6), "All windows fell have coordinates outside chromosome boundaries")
		
	})

# ---------------------------------------------------------------------------- #
test_that("scoreMatrix:GRanges, GRanges works",
{
          target = GRanges(rep(c(1,2),each=6), IRanges(rep(c(1,2,3,7,8,9), times=2), width=5))
          values(target)$weight = rep(c(1,2),each=6)
          windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5))
          
          
          # normal function works
          s1 = ScoreMatrix(target=target, windows=windows)
          m1 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
          rownames(m1) = 1:4
          m1 = as(m1, 'ScoreMatrix')
          expect_equal(s1,m1)
          
          # function with weight col works
          s2 = ScoreMatrix(target=target, windows=windows, weight.col='weight')
          m2 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
          m2[3:4,] = m2[3:4,]*2
          rownames(m2) = 1:4
          m2 = as(m2, 'ScoreMatrix')
          expect_equal(s2,m2)      
})


# ---------------------------------------------------------------------------- #
test_that("scoreMatrix:character, GRanges works",
{
  target = GRanges(rep(c(1,2),each=6), IRanges(rep(c(1,2,3,7,8,9), times=2), width=5))
  values(target)$weight = rep(c(1,2),each=6)
  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5))
  
 
  
  # function with weight col works
  s2 = ScoreMatrix(target=target, windows=windows, weight.col='weight')
  m2 = matrix(rep(c(1,2,3,3,3,2,3,3,3,2),times=2), ncol=5, byrow=T)
  m2[3:4,] = m2[3:4,]*2
  rownames(m2) = 1:4
  m2 = as(m2, 'ScoreMatrix')
  expect_equal(s2,m2)
  
  
  
})
	
# ---------------------------------------------------------------------------- #
# test for removeOffRanges
test_that("constrainRanges work",
	{
		target= RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,8),c(3,7,3,10)))
		gt1 = GRanges(c('chr1','chr1','chr2'), IRanges(c(1,5,1),c(3,7,3)))
		expect_identical(constrainRanges(target, gr1), gt1)
	}
)

# ---------------------------------------------------------------------------- #
# test for binMatrix
test_that("binMatrix works"
	m1 = as(matrix(rep(1:4, 4), ncol=4), 'scoreMatrix')
	
	# no binning
	expect_identical(binMatrix(m1), m1)
	
	# nbins 2, default function
	mb1 = as(matrix(c(1,1,2,2,3,3,4,4), ncol=2, byrow=T), 'scoreMatrix')
	expect_identical(binMatrix(m1, nbins=2), mb1)
	
	# nbins 5
	expect_error(binMatrix(m1, nbins=5), "number of given bins is bigger than the number of matrix columns")
	
	# nbins 2, not default function
	m2 = as(matrix(rep(1:4, 4), ncol=4, byrow=T), 'scoreMatrix')
	mb2 = as(matrix(rep(c(2,4), 4), ncol=2, byrow=T), 'scoreMatrix')
	expect_equal(binMatrix(m2, nbins=2, max), mb2)	
)

# ---------------------------------------------------------------------------- #
test_that("ScoreMatrix BigWig works"
	library(rtracklayer)
	test_path <- system.file("tests", package = "rtracklayer")
	test_bw <- file.path(test_path, "test.bw")
	b = import(test_bw, asRangedData=F)

	st = seq(200, 300, 20)
	g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
	s = ScoreMatrix(test_bw, g, type='bigWig')

	m = matrix(-1, ncol=10, nrow=6)
	m[6,-1] = -0.75
	rownames(m) = 1:6
	m = as(m, 'ScoreMatrix')
	expect_equal(s, m)
)


# ---------------------------------------------------------------------------- #
test_that("scaleScoreMatrix works"
          s = as(matrix(1:9, ncol=3),'ScoreMatrix')
          
          s1=scaleScoreMatrix(s, rows=T, columns=F)
          s2=scaleScoreMatrix(s, rows=F, columns=T)
          s3=scaleScoreMatrix(s, rows=T, columns=T)
          
          m1 = as(matrix(c(rep(-0.4285714,3),rep(0,3), rep(0.4285714,3)), 
                         nrow=3),'ScoreMatrix')
          m2 = as(matrix(c(rep(-0.3333333,3),rep(0,3), rep(0.3333333,3)), 
                         nrow=3, byrow=T),'ScoreMatrix')
          m3 = as(matrix(0,nrow=3,ncol=3),'ScoreMatrix')
          
          expect_equal(s1,m1)
          expect_equal(s2,m2)
          expect_equal(s3,m3)
)


        