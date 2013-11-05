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
		
		#7. for modRle
		gm = GRanges(rep(c('chr1','chr2'), each=4), IRanges(rep(2:5, times=2), rep(3:6, times=2)))
		values(gm)$weight = rep(2:3, 4)
		mc = modCoverage(gm, 'weight', multiply=1, add=0)
		grm = GRanges(rep(c('chr1','chr2'), each=2), IRanges(rep(c(1,3), 2), rep(c(3,5), 2)))
		
		m7 = as(matrix(c(0, 2, 5, 5, 5, 5, 0, 2, 5, 5, 5, 5), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m7) = 1:4
		expect_equal(scoreMatrix(mc, grm), m7)
	})
	
# ---------------------------------------------------------------------------- #
# test for removeOffRanges
test_that("constrainRanges work",
	{
		target= RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,8),c(3,7,3,10)))
		gt1 = GRanges(c('chr1','chr1','chr2'), IRanges(c(1,5,1),c(3,7,3)))
		expect_identical(removeOffRanges(target, gr1), gt1)
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

	s = seq(1, 1000, 20)
	g = GRanges(rep('chr2', length(s)), IRanges(s, width=50))
	s = ScoreMatrix(test_bw, g, type='bigWig')

	m = matrix(-1, ncol=10, nrow=6)
	m[6,-1] = -0.75
	rownames(m) = 1:6
	m = as(m, 'ScoreMatrix')
	expect_equal(s, m)
)
