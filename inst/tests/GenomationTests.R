# ------------------------------------------- #
# test for removeOffRanges
test_that("constrainRanges work",
	{
		target= RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,8),c(3,7,3,10)))
		gt1 = GRanges(c('chr1','chr1','chr2'), IRanges(c(1,5,1),c(3,7,3)))
		expect_identical(removeOffRanges(target, gr1), gt1)
	}
)



# ------------------------------------------- #
# test for scoreMatrix function
test_that("scoreMatrix works",
	{
	
		# input RleList
		rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		
		#1. test for proper workings
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)), strand=c('+','-','+','-'))
		m1 = as(matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m1) = 1:4
		expect_identical(scoreMatrix(rl, gr1, strand.aware = FALSE), m1)
		
		m2 = as(matrix(c(1,1,1,3,2,2,4,4,4,6,5,5), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m2) = 1:4
		expect_identical(scoreMatrix(rl, gr1, strand.aware = TRUE), m2)
		
		#2. test for different lengths
		gr2 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,9)))
		expect_error(scoreMatrix(rl, gr2), "width of 'windows' are not equal, provide 'windows' with equal widths")
		
		#3. test for removing outliers
		gr3 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,8,1),c(3,7,10,3)))
		m3 = as(matrix(c(1,1,1,2,2,3,4,4,4), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m3) = c(1,2,4)
		expect_that(scoreMatrix(rl, gr3), is_identical_to(m3))
		
		#4. test for different chromosomes
		gr4 = GRanges(rep(c('chr1','chr3'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		m4 = as(matrix(c(1,1,1,2,2,3), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m4) = 1:2
		expect_identical(scoreMatrix(rl, gr4, strand.aware = FALSE), m4)
		
		#5. no overlapping chromosomes
		gr5 = GRanges(rep(c('chr3','chr4'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		expect_error(scoreMatrix(rl, gr5), "all windows fell have coordinates outside chromosome boundaries")
		
		#6. only windows that fall off chromosomes
		gr6 = GRanges(rep(c('chr1','chr4'), each=2), IRanges(c(-1,8,1,5),c(1,10,3,7)))
		expect_warning(scoreMatrix(rl, gr6), "all windows fell have coordinates outside chromosome boundaries")
		
		#7. for modRle
		gm = GRanges(rep(c('chr1','chr2'), each=4), IRanges(rep(2:5, times=2), rep(3:6, times=2)))
		values(gm)$weight = rep(2:3, 4)
		mc = modCoverage(gm, 'weight', multiply=1, add=0)
		grm = GRanges(rep(c('chr1','chr2'), each=2), IRanges(rep(c(1,3), 2), rep(c(3,5), 2)))
		
		m7 = as(matrix(c(0, 2, 5, 5, 5, 5, 0, 2, 5, 5, 5, 5), ncol=3, byrow=T), 'scoreMatrix')
		rownames(m7) = 1:4
		expect_that(scoreMatrix(mc, grm), is_identical_to(m7))	
	})
	

# ------------------------------------------- #
# test for scoreMatrixBin function
test_that("scoreMatrixBin works",
	{
	
		# input RleList
		rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		
		#1. test for proper workings
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(2,7,4,9)), strand=c('+','-','+','-'))
		m1 = as(matrix(c(1,1,2,2.5,4,4.5,5.3333333333,6), ncol=2, byrow=T), 'scoreMatrix')
		rownames(m1) = 1:4
		expect_that(scoreMatrixBin(rl, gr1, strand.aware = FALSE, bin.num=2), equals(m1))
		
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
	

# ------------------------------------------- #
# test for binner function
test_that("binner works",
	{
		#1.
		m1 = matrix(c(1,2))
		rownames(m1) = c('my.start','my.end')
		expect_that(binner(1,2,1), is_equivalent_to(m1))
		
		#2.
		m2 = matrix(c(1,1,2,2), ncol=2)
		rownames(m2) = c('my.start','my.end')
		expect_that(binner(1,2,2), is_equivalent_to(m2))
		
		#3.
		m3 = matrix(c(1,1,2,1,2,2), ncol=3)
		rownames(m3) = c('my.start','my.end')
		expect_that(binner(1,2,3), is_equivalent_to(m3))	
	})
	
	