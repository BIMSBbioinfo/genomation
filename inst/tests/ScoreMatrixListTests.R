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
test_that("ScoreMatrixList BigWig works"
	library(rtracklayer)
	test_path <- system.file("tests", package = "rtracklayer")
	test_bw <- file.path(test_path, "test.bw")

  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
	s = ScoreMatrixList(c(test_bw, test_bw), g, type='bigWig')

	m = matrix(-1, ncol=10, nrow=6)
	m[6,-1] = -0.75
	rownames(m) = 1:6
	m = as(m, 'ScoreMatrix')
	l = as(list(test.bw=m, test.bw=m), 'ScoreMatrixList')
	expect_equal(s, l)
)

test_that("ScoreMatrixList Bin BigWig works"
	library(rtracklayer)
	test_path <- system.file("tests", package = "rtracklayer")
	test_bw <- file.path(test_path, "test.bw")

  st = seq(200, 300, 20)
  g = GRanges(rep('chr2', length(st)), IRanges(st, width=10))
  b = import(test_bw, asRangedData=FALSE, which=g)
  covs = coverage(b, weight=b$score)        
  s = ScoreMatrixBin(covs, g, bin.num=5)
          
  s = ScoreMatrixList(c(test_bw, test_bw), g,bin.num=5, type='bigWig')

	m = matrix(-1, ncol=10, nrow=6)
	m[6,-1] = -0.75
	rownames(m) = 1:6
	m = as(m, 'ScoreMatrix')
  m = binMatrix(m, bin.num=5)
	l = as(list(test.bw=m, test.bw=m), 'ScoreMatrixList')
	expect_equal(s, l)
)
