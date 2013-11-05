# ---------------------------------------------------------------------------- #
# test for scoreMatrixBin function
test_that("scoreMatrixBin works",
	{
	
		# input RleList
		rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		
		#1. test for proper workings
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(2,7,4,9)), strand=c('+','-','+','-'))
		m1 = as(matrix(c(1,1,2,2.5,4,4.5,5.3333333333,6), ncol=2, byrow=T), 'scoreMatrix')
		rownames(m1) = 1:4
		expect_equal(scoreMatrixBin(rl, gr1, strand.aware = FALSE, bin.num=2), m1)
		
		#2. test for windows smaller than bin.num
		m2 = as(matrix(c(2,2,3.0,4,4,4.5, 5,6.0,6), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m2) = 2:4
		expect_that(scoreMatrixBin(rl, gr1, strand.aware = FALSE, bin.num=3), equals(m2))
		expect_that(scoreMatrixBin(rl, gr1, strand.aware = FALSE, bin.num=3), gives_warning("supplied GRanges object contains ranges of width < number of bins"))

		#3. for modRle
		gm = GRanges(rep(c('chr1','chr2'), each=8), IRanges(rep(2:9, times=2), rep(3:10, times=2)))
		values(gm)$weight = rep(2:3, 8)
		mc = modCoverage(gm, 'weight', multiply=1, add=0)
		grm = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(2,7,4,9)), strand=c('+','-','+','-'))
		# scoreMatrixBin(mc, grm, 2)	
	})