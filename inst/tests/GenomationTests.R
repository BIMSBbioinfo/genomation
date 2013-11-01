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
# scoreMatrixList tests
l = lapply(seq(20, 40,5), function(x)matrix(rpois(1000, x), ncol=25))
l = lapply(l, function(x)as(x, "scoreMatrix"))
l = scoreMatrixList(l)
png(file.path(getwd(),'sml.png'), width=800, height=800)
	heatmapProfile(l, xcex=1.5, ycex=1.5, cex.main=3,  ylab='Number',  xlab='Position')
dev.off()
	
# tests the yargs and y.at arguments
load_all(pkg = genomation.path)
l = lapply(seq(20, 40,5), function(x)matrix(rpois(1000, x), ncol=25))
l = lapply(l, function(x)as(x, "scoreMatrix"))
l = scoreMatrixList(l)
png(file.path('/home/members/vfranke/Tmp','sml1.png'), width=800, height=800)
	ymarks1 = c('aa','bb','cc','dd')
	y.at1 = c(10,20,30,40)	
	heatmapProfile(l, xcex=1.5, ycex=1.5, cex.main=3,  ylab='Number',  xlab='Position', ymarks=ymarks1, y.at=y.at1)
dev.off()
	
png(file.path('/home/members/vfranke/Tmp','sml2.png'), width=1200, height=800)
	ymarks2 = c('aaaaaaaaaaaaaaaaaaaaaaaaa','bb','cc','dd')
	y.at2 = c(10,12,14,16)	
	heatmapProfile(l, xcex=1.5, ycex=1.5, cex.main=3,  ylab='Number',  xlab='Position', ymarks=ymarks2, y.at=y.at2)
dev.off()

# ---------------------------------------------------------------------------- #
# test for plotMatrix
m = as(matrix(rnorm(500), nrow=50), 'scoreMatrix')
rownames(m) = sample(letters, nrow(m), replace=T)
heatMatrix(m, use.names=T)
png(file.path('/home/members/vfranke/Tmp','pm.png'), width=1200, height=1200)
	
dev.off()


	
	