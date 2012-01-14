### INFO: Development tools for the genomation package
### DATE: 08.11.2011.
### AUTHOR: frenkiboy

# {1} LIBRARIES
library(devtools)
 library(testthat)
#/{1} LIBRARIES


# {2} CODE
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		lib.path='/home/members/vfranke/Projects/Code/Scripts/TestLib'
		genomation.path='/home/members/vfranke/Projects/Code/Scripts/Genomation'
		
		
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	dev_mode(on=TRUE, path=lib.path)
	install(genomation.path)
	load_all(pkg = genomation.path)

	code_path = file.path(genomation.path, 'R')
	test_path = file.path(genomation.path, 'inst', 'tests')
	test_dir(test_path)
	
	set.seed(10)
	
	f = function(x, n=2000){
		n=n
		r = RleList('1'=c(sapply(c(x, x+5, x+10, x+15), function(y)rpois(n, y))))
		g = GRanges(1, IRanges(seq(1, length(r[[1]]), 50), width=50))
		scoreMatrix(r, g)
	}
	l = lapply(c(5,10,15,20), f)
	lg = scoreMatrixList(l)
	
	png(file.path(lib.path, 'sml.png'), width=1000, height=800)
		heatmapProfile(lg)
	dev.off()
		
	
	sl1 = s1[1:4,1:4]
	sl2 = s1[1:4,1:10]
	sl3 = cbind(sl2,sl2)
	sl4 = rbind(sl3, sl3, sl3)
	sl = scoreMatrixList(sl1,sl2,sl3,sl4)
	s2 = scoreMatrix(r, g)
	sl = scoreMatrixList(s1, s2)
	
	#/{{3}} MAIN
#/{2} CODE



