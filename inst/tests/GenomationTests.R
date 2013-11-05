
# ---------------------------------------------------------------------------- #
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
	



# ---------------------------------------------------------------------------- #
# test for plotMatrix
m = as(rbind(matrix(rnorm(250), nrow=50),matrix(rnorm(250, 5), nrow=50)), 'ScoreMatrix')
rownames(m) = sample(letters, nrow(m), replace=T)
heatMatrix(m, use.names=T)
heatMatrix(m, use.names=F)
fact = factor(rep(c(1,2), each=50))
heatMatrix(m, use.names=F, fact=fact)


	
	