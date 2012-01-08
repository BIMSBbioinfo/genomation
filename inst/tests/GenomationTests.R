# test for scoreMatrix function
test_that("scoreMatrix works",
	{
	
		# input RleList
		rl = RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		
		# test for proper workings
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,7)))
		m1 = as(matrix(c(1,1,1,2,2,3,4,4,4,5,5,6), ncol=3, byrow=T), 'scoreMatrix')
		expect_identical(scoreMatrix(rl, gr1), m1)
		
		# test for different lengths
		gr2 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,5),c(3,7,3,9)))
		expect_error(scoreMatrix(rl, gr2), "width of 'windows' are not equal, provide 'windows' with equal widths")
		
		# test for removing outliers
		gr3 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,8),c(3,7,3,10)))
		m3 = as(matrix(c(1,1,1,2,2,3,4,4,4), ncol=3, byrow=T), 'scoreMatrix')
		expect_that(scoreMatrix(rl, gr3), is_identical_to(m3))
	})
