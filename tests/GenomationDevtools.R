### INFO: Development tools for the genomation package
### DATE: 08.11.2011.
### AUTHOR: frenkiboy

# {1} LIBRARIES
library(devtools)
#/{1} LIBRARIES


# {2} CODE
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		lib.path='/home/members/vfranke/Projects/Code/Scripts/TestLib'
		genomation.path='/home/members/vfranke/Projects/Code/Scripts/Genomation/Genomation'
		
		
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	dev_mode(on=TRUE, path=lib.path)
	install(genomation.path)
	load_all(pkg = genomation.path)
	
	#/{{3}} MAIN
#/{2} CODE



