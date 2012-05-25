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
	peaks = lapply(peak.files, read.table, header=F)
	names(peaks) = str_replace(basename(peak.files),'.narrowPeak','')
	peaks = lapply(peaks, function(x)convert.bed.df(x[,1:3]))

	genes = read.transcript.features(genes.path)
	annot = lapply(peaks, function(x)annotate.WithGenicParts(x, genes))
	
	library(ggplot2)
	library(reshape2)
	d = do.call(rbind, lapply(annot, function(x)x@precedence))
	h = hclust(dist(d))
	# d = cbind(d, id=1:nrow(d), ids=order((1:nrow(d))[h$order]))
	m = reshape2::melt(data.frame(d))
	m = melt(d)
	
	# CairoPNG(file.path(lib.path, "heatmap.png"), width=800, height=800)
		# pheatmap(mat=d, kmeans_k=NA, scale="none", cluster_rows=T, cluster_cols=F, clustering_distance_rows='euclidean')
	# dev.off()
	
	CairoPNG(file.path(lib.path, "heatmap.png"), width=800, height=800)
		p <- ggplot(m, aes(x=X2, y=X1, fill=value, colour="white")) + scale_fill_gradient(low = "white",high = "cornflowerblue") + scale_y_discrete(limits = rownames(d)[h$order] ) + opts(axis.title.x=theme_text(colour='white'), axis.text.x=theme_text(colour='black', face='bold'), axis.text.y=theme_text(colour='black'), axis.title.y=theme_text(colour='white', face='bold'))
		p + geom_tile(color='white') 
	dev.off()
	
	library(Cairo)
	library(RColorBrewer)
	pal = brewer.pal(9, 'Blues')
	CairoPNG(file.path(lib.path, "heatmap.png"), width=800, height=800)
		p <- ggplot(m, aes(x=X2, y=X1, fill=value, colour="white")) + scale_fill_gradient(low = "white",high = "steelblue") + scale_x_discrete() + opts(axis.title.x=theme_text(colour='white'), axis.text.x=theme_text(colour='black', face='bold'), axis.text.y=theme_text(colour='black'), axis.title.y=theme_text(colour='white', face='bold'))
		p + geom_tile(color='white') 
	dev.off()
	+ scale_fill_manual(values=my.cols) + scale_x_discrete(labels=lab.lookup$feature.labels)
	
		# opts(axis.text.x=theme_text(hjust=1, angle=90),
       # axis.text.y=theme_text(hjust=1)) + scale_colour_discrete(guide="none") + opts(legend.position = "top")
	
	
	#/{{3}} MAIN
#/{2} CODE



