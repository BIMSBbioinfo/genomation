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
	n = 2000
	r = RleList('1' = (c(sapply(c(5,15, 25, 35), function(x)rpois(n, x)))))
	g = GRanges(1, IRanges(seq(1, length(r[[1]]), 50), width=50))
	s = scoreMatrix(r, g)
	f = factor(rep(1:4, each=n/50))
	png(file.path(lib.path, 'scoreMatrix.png'), width=1000, height=600)
		plotMatrix(s, fact=f)
	dev.off()
	png(file.path(lib.path, 'scoreMatrix.png'), width=1000, height=600)
		plotMatrix(s)
	dev.off()
	
	levels(f) = c('Class1','Class2','Class3','Class4')
	png(file.path(lib.path, 'scoreMatrix.png'), width=1000, height=600)
		plotMatrix(s, fact=f)
	dev.off()
	
	#/{{3}} MAIN
#/{2} CODE



