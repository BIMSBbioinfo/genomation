### INFO: Development tools for the genomation package
### DATE: 08.11.2011.
### AUTHOR: frenkiboy

# {1} LIBRARIES
library(devtools)
 # library(testthat)
#/{1} LIBRARIES


# {2} CODE
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		lib.path='/home/members/vfranke/Projects/Code/Genomation/TestLib'
		genomation.path='/home/members/vfranke/Projects/Code/Genomation/genomation'
		
		
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
	
	genes.path = '/common/USERS/vfranke/Work/Genomation/Data/Ensembl/hg19.ensembl.bed'
	peaks.path = '/common/USERS/vfranke/Work/Genomation/Data/Encode/SYDH/Gm12878/Peaks'
	rdata.path = ''

	peak.files = list.files(peaks.path, pattern='narrow', full.names=T)
	peaks = lapply(peak.files, read.table, header=T)
	names(peaks) = str_replace(basename(peak.files),'.narrowPeak','')

	genes = read.transcript.features(genes.path)
	
	
	#/{{3}} MAIN
#/{2} CODE



