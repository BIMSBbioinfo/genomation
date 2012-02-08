# ------------------------------------------- #
# test for removeOffRanges
test_that("removeOffRanges workd",
	{
		target= RleList(chr1 = Rle(rep(c(1,2,3), each=3)), chr2=Rle(rep(c(4,5,6), each=3)))
		gr1 = GRanges(rep(c('chr1','chr2'), each=2), IRanges(c(1,5,1,8),c(3,7,3,10)))
		gt1 = GRanges(c('chr1','chr1','chr2'), IRanges(c(1,5,1),c(3,7,3)))
		expect_identical(removeOffRanges(target, gr1), gt1)
	}
)



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
		expect_error(scoreMatrix(rl, gr5), "There are no common chromosomes/spaces to do overlap")
		
		#6. only windows that fall off chromosomes
		gr6 = GRanges(rep(c('chr1','chr4'), each=2), IRanges(c(-1,8,1,5),c(1,10,3,7)))
		expect_warning(scoreMatrix(rl, gr6), "windows have no ranges left after filtering")
		
	})
	

# ------------------------------------------- #
# test for scoreMatrixBin function
test_that("scoreMatrixBin works",
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
		expect_error(scoreMatrix(rl, gr5), "There are no common chromosomes/spaces to do overlap")
		
		#6. only windows that fall off chromosomes
		gr6 = GRanges(rep(c('chr1','chr4'), each=2), IRanges(c(-1,8,1,5),c(1,10,3,7)))
		expect_warning(scoreMatrix(rl, gr6), "windows have no ranges left after filtering")
		
	})
	

# TestScaler = function(){
  
  # error = character()
    
  # test1 
  # r1 = Rle(c(1,2,3,4,5,5,5,3,2))
  # l1 = 2
  # if(!all(ScalerLarge(r1, l1) == c(3,4))){
    # error = c(error, paste('Error:', length(r1), l1, "\n"))
  # }
    
  # test2
  # r2 = Rle(c(1,2,3,5,5,5,4,4,4))
  # l2 = 3
  # if(!all(ScalerLarge(r2, l2) == c(2,5,4))){
    # error = c(error, paste('Error:', length(r2), l2, "\n"))
  # }

  # r3 = Rle(c(2,2,2,3,4,5,6,7,8))
  # l3 = 4
  # if(!all(ScalerLarge(r3, l3) == c(2,3,5,7))){
    # error = c(error, paste('Error:', length(r3), l3, "\n"))
  # }
 
  # r4 = Rle(c(2,2,3,3,4,5,5,6,6))
  # l4 = 5
  # if(!all(ScalerLarge(r4, l4) == c(2,3,4,5,6))){
    # error = c(error, paste('Error:', length(r4), l4, "\n"))
  # }
 
  # if(length(error) == 0){
    # return(TRUE)
  # }else{
    # stop(error)
  # } 
# }

	
# ------------------------------------------- #
# test for scoreMatrixList
# test_that("scoreMatrixList initalization works",
	# {
		
	
	# })
	