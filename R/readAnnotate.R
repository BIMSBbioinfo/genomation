# PART THAT DEALS with annotation of genomic features

##############################################################################
# SECTION 1:
# reading annotation to GRanges
# makes GRanges object from a given bed6 or bed12 file to granges object
##############################################################################

#######################################
# SECTION 1: S3 functions
#######################################

# ----------------------------------------------------------------------------------------------- #
# extracts exons from a bed12 file and puts them into GRanges object
bed12ToExons<-function(ref){

	ref=unique(ref)

	# apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
	b.start.size = cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),as.integer(unlist(strsplit(as.character(ref$V11),"," ))))
	rep.ref = ref[rep(1:nrow(ref),ref[,10]),] # replicate rows occurs as many as its exon number


	exon.id = unlist( mapply( function(x,y) if(x=="+"){return(1:y)}else{return(y:1)} ,ref[,6],ref[,10]  ) )
	rep.ref$V5=exon.id

	rep.ref$V3 = rep.ref$V2+b.start.size[,1]+b.start.size[,2] # calculate exon start and ends
	rep.ref$V2 = rep.ref$V2+b.start.size[,1]

	strands = as.character(rep.ref$V6)
	strands[strands=="."]="*"
	grange.exons = GenomicRanges::GRanges(seqnames=as.character(rep.ref$V1),
	ranges = IRanges(start=rep.ref$V2+1, end=rep.ref$V3),
	strand = strands, 
	score=rep.ref$V5,
	name=rep.ref$V4)
	return(grange.exons)
} 

# ----------------------------------------------------------------------------------------------- #
# extracts introns from a bed12 file and puts them into GRanges object
bed12ToIntrons<-function(ref){

	#remove the genes with one exon only (they won't have any introns)
	ref=ref[ref[,10]>1,]
	ids=paste(ref[,1],ref[,2],ref[,3],ref[,4],sep="")
	ref=cbind(ref,id=ids)
	ref=unique(ref)

	# apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
	b.start.size=cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),as.integer(unlist(strsplit(as.character(ref$V11),"," ))))

	# replicate rows occurs as many as its exon number
	rep.ref = ref[rep(1:nrow(ref),ref[,10]),] 

	exon.id = unlist( mapply( function(x,y) if(x=="+"){return(1:y)}else{return(y:1)} ,ref[,6],ref[,10]  ) )
	rep.ref$V5 = exon.id

	# calculate exon start and ends
	rep.ref$V3 = rep.ref$V2+b.start.size[,1]+b.start.size[,2] 
	rep.ref$V2 = rep.ref$V2+b.start.size[,1]
	rep.ref = rep.ref[,c(1:6,13)]

	# now try to get the exons by cbinding consecutive exons
	temp.ref		= cbind(rep.ref[1:(nrow(rep.ref)-1),],rep.ref[2:nrow(rep.ref),])
	temp.ref		= temp.ref[temp.ref[,7]==temp.ref[,14],]
	temp.ref[,2]	= temp.ref[,3]
	temp.ref[,3]	= temp.ref[,9]
	rep.ref			= temp.ref[,1:6]

	# subtract 1 from - strand exon numbers
	rep.ref[rep.ref[,6]=="-",5]=rep.ref[rep.ref[,6]=="-",5]-1 

	strands=as.character(rep.ref$V6)
	strands[strands=="."]="*"
	grange.exons=GenomicRanges::GRanges(seqnames=as.character(rep.ref$V1),
	ranges=IRanges(start=rep.ref$V2+1, end=rep.ref$V3),
	strand=strands, score=rep.ref$V5,name=rep.ref$V4)
	return(grange.exons)
} 

# ----------------------------------------------------------------------------------------------- #
# checks the validity of the bed data.frame if it is a legitimate bed columns
checkBedValidity<-function(bed.df,type="none"){

	# does it have at least 3 columns
	num.col=(ncol(bed.df)>=3)

	# does it have 1st and 2nd column as numbers
	col1.2= (is.integer(bed.df[,2]) & is.integer(bed.df[,3]) )

	# does it have chr string in 1st column
	chr= sum(grepl("chr",bed.df[,1]) )

	if(type=="exons" | type=="introns") {
		#does it have 12>= columns
		ex=(ncol(bed.df)>=12) 

		return(num.col & col1.2 & chr & ex)
	}
	return(num.col & col1.2 & chr )
}




#######################################
# SECTION 1: S4 functions
#######################################


# ----------------------------------------------------------------------------------------------- #
#' convert a data frame read-in from a bed file to a GRanges object
#'  
#' @param bed  a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convertBedDf(bed)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#' @note bed files are expected not to have header lines
#'
#' @export
#' @docType methods
#' @rdname convertBedDf-methods
setGeneric("convertBedDf",function(bed) standardGeneric("convertBedDf"))

#' @aliases convertBedDf,data.frame-method
#' @rdname convertBedDf-methods
setMethod("convertBedDf" ,
		  signature(bed = "data.frame" ),
		  function(bed){

			if(! checkBedValidity(bed))
				stop("this is not a valid bed file")
				
			if(ncol(bed)>=6){
			
				# check if 6th column is really strand
				if( sum( unique(bed[,6]) %in% c("+","-",".") ) != length(unique(bed[,6])) )
					stop("Strand column of the bed file or data.frame is wrong")
						
				#convert to granges
				strands=as.character(bed$V6)
				strands[strands=="."]="*"
				grange=GRanges(seqnames=as.character(bed$V1),
							   ranges=IRanges(start=bed$V2+1, end=bed$V3),
							   strand=strands, 
							   score=bed$V5,
							   name=bed$V4)
			}

			if(ncol(bed)==5){
				grange = GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),strand=rep("*",nrow(bed)), score=bed$V5,name=bed$V4 )
			}
			if(ncol(bed)==4){
				grange = GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),strand=rep("*",nrow(bed)),name=bed$V4)
			}    
			if(ncol(bed)==3){
				grange = GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),strand=rep("*",nrow(bed)))
			}                           
			return(grange)
})

# ----------------------------------------------------------------------------------------------- #
#' convert a data frame read-in from a bed file to a GRanges object for exons
#'  
#' @param bed.df  a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convertBed2Exons(bed.df)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname convertBed2Exons-methods
setGeneric("convertBed2Exons",function(bed.df) standardGeneric("convertBed2Exons"))

#' @aliases convertBed2Exons,data.frame-method
#' @rdname convertBed2Exons-methods
setMethod("convertBed2Exons" ,
			signature(bed.df = "data.frame" ),
			function(bed.df){

			if(! checkBedValidity(bed.df,"exon"))
				stop("this is not a valid bed file")
				
			bed12ToExons(bed.df)
})


# ----------------------------------------------------------------------------------------------- #
#' convert a data frame read-in from a bed file to a GRanges object for introns
#'  
#' @param bed.df  a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convertBed2Introns(bed.df)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname convertBed2Introns-methods
setGeneric("convertBed2Introns",function(bed.df) standardGeneric("convertBed2Introns"))

#' @aliases convertBed2Introns,data.frame-method
#' @rdname convertBed2Introns-methods
setMethod("convertBed2Introns",
			signature(bed.df = "data.frame" ),
			function(bed.df){

			if(! checkBedValidity(bed.df,"exon"))
				stop("this is not a valid bed file")
				
			bed12ToIntrons(bed.df)
})

# ----------------------------------------------------------------------------------------------- #
#' Function for reading exon intron and promoter structure from a given bed file
#'
#' @param location location of the bed file with 12 or more columns
#' @param remove.unsual remove the chromomesomes with unsual names, mainly random chromsomes etc
#' @param up.flank  up-stream from TSS to detect promoter boundaries
#' @param down.flank down-stream from TSS to detect promoter boundaries
#' @param unique.prom     get only the unique promoters, promoter boundaries will not have a gene name if you set this option to be TRUE
#' @usage readTranscriptFeatures(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE)
#' @return a \code{\link{GRangesList}} containing locations of exon/intron/promoter/TSS
#' @note  one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#' 
#' @examples
#'   my.bed12.file=my.file=system.file("extdata", "hg18.refseq.txt.Test", package = "genomation")
#'   my.bed12.file
#'   feats=readTranscriptFeatures(my.bed12.file) 
#'
#' @export
#' @docType methods
#' @rdname readTranscriptFeatures-methods
setGeneric("readTranscriptFeatures", function(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE) standardGeneric("readTranscriptFeatures"))

#' @aliases readTranscriptFeatures,character-method
#' @rdname readTranscriptFeatures-methods
setMethod("readTranscriptFeatures", 
			signature(location = "character"),#,remove.unsual="logical",up.flank="numeric",down.flank="numeric",unique.prom="logical" ),
			function(location,remove.unsual,up.flank ,down.flank ,unique.prom){

			# find out if there is a header, skip 1st line if there is a header
			f.line=readLines(con = location, n = 1)
			skip=0
			if(grepl("^track",f.line))
				skip=1

			# readBed6
			cat('Reading the table...\r')
			bed=.readTableFast(location,header=F,skip=skip)                    
			if(remove.unsual)
				bed=bed[grep("_", as.character(bed[,1]),invert=T),]
			
			# introns
			cat('Calculating intron coordinates...\r')
			introns	= convertBed2Introns(bed)
			# exons
			cat('Calculating exon coordinates...\r')
			exons	= convertBed2Exons(bed)

			# get the locations of TSSes
			cat('Calculating TSS coordinates...\r')
			tss=bed
			#  + strand
			tss[tss$V6=="+",3] = tss[tss$V6=="+",2]
			#  - strand
			tss[tss$V6=="-",2]=tss[tss$V6=="-",3]

			tssg = GRanges(seqnames=as.character(tss$V1),
						   ranges=IRanges(start=tss$V2, end=tss$V3),
						   strand=as.character(tss$V6),
						   score=rep(0,nrow(tss)),
						   name=tss$V4)

			cat('Calculating promoter coordinates...\r')
			# get the locations of promoters
			# + strand
			bed[bed$V6=="+",3]=bed[bed$V6=="+",2]+down.flank
			bed[bed$V6=="+",2]=bed[bed$V6=="+",2]-up.flank

			#  - strand
			bed[bed$V6=="-",2]=bed[bed$V6=="-",3]-down.flank
			bed[bed$V6=="-",3]=bed[bed$V6=="-",3]+up.flank

			
			if(! unique.prom){
				prom.df = (bed[,c(1,2,3,4,6)])
				prom = GRanges(seqnames=as.character(prom.df$V1),
							   ranges = IRanges(start=prom.df$V2, end=prom.df$V3),
							   strand = as.character(prom.df$V6),
							   score=rep(0,nrow(prom.df)),
							   name=prom.df$V4)
			}else{
				prom.df = unique(bed[,c(1,2,3,6)])
				prom = GRanges(seqnames=as.character(prom.df$V1),
							   ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
							   strand=as.character(prom.df$V6),
							   score=rep(0,nrow(prom.df)),
							   name=rep(".",nrow(prom.df)) )
			}
			
			cat('Outputting the final GRangesList...\r\n')
			GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
})


# ----------------------------------------------------------------------------------------------- #
#' Function to get upstream and downstream adjecent regions to a genomic feature such as CpG islands
#' 
#' @param grange GRanges object for the feature
#' @param flank  number of basepairs for the flanking regions
#' @param clean  If set to TRUE, flanks overlapping with other main features will be trimmed, and overlapping flanks will be removed
#'        this will remove multiple counts when other features overlap with flanks
#'
#' @usage getFlanks(grange,flank=2000,clean=T)
#' @return GRanges object for flanking regions
#' @export
#' @docType methods
#' @rdname getFlanks-methods
setGeneric("getFlanks", function(grange,flank=2000,clean=T) standardGeneric("getFlanks"))

#' @aliases getFlanks,GRanges-method
#' @rdname getFlanks-methods
setMethod("getFlanks", signature(grange= "GRanges"),
			function(grange,flank=2000,clean=T){

			shores = c( IRanges::flank(grange,flank),
						IRanges::flank(grange,flank, FALSE) )
			
			# merge overlapping shores remove CpG coordinates from all shores, das ist so cool!!
			if(clean)
				shores=IRanges::reduce(IRanges::setdiff(shores, grange)) 
			
			shores
})


# ----------------------------------------------------------------------------------------------- #
#' A function to read-in genomic features and their upstream and downstream adjecent regions such as CpG islands and their shores
#'
#' @param location for the bed file of the feature 
#' @param flank    number of basepairs for the flanking regions
#' @param clean    If set to TRUE, flanks overlapping with other main features will be trimmed
#' @param remove.unsual  remove chromsomes with unsual names random, Un and antyhing with "_" character
#' @param feature.flank.name the names for feature and flank ranges, it should be a character vector of length 2. example: c("CpGi","shores")
#' @usage  readFeatureFlank(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL)
#' @return a GenomicRangesList contatining one GRanges object for flanks and one for GRanges object for the main feature.
#'   NOTE:This can not return a GRangesList at the moment because flanking regions do not have to have the same column name as the feature.
#'   GRangesList elements should resemble eachother in the column content. We can not satisfy that criteria for the flanks
#'
#' @export
#' @docType methods
#' @rdname readFeatureFlank-methods
setGeneric("readFeatureFlank", function(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL) standardGeneric("readFeatureFlank") )

#' @aliases readFeatureFlank,character-method
#' @rdname readFeatureFlank-methods
setMethod("readFeatureFlank", 
			signature(location = "character"),
			function(location,remove.unsual,flank ,clean,feature.flank.name){

				feat = readBed(location, remove.unsual)
				flanks = getFlanks(feat,flank=flank,clean=clean)
				x = GenomicRangesList(features=feat,flanks=flanks)
				if(!is.null(feature.flank.name) & length(feature.flank.name)==2)
					names(x)=feature.flank.name

				return(x)
})



##############################################################################
# SECTION 2:
# annotate granges objects with annotations that read-in and converted to GRanges objects
##############################################################################

#######################################
# SECTION 2: Define new classes
#######################################

# ----------------------------------------------------------------------------------------------- #
# A set of objects that will hold statistics about feature and annotation overlap
#' An S4 class that information on overlap of target features with annotation features  
#'
#' This object is desgined to hold statistics and information about genomic feature overlaps
#'          
#' @section Slots:\describe{
#'            	    \item{\code{members}}{a matrix showing overlap of target features with annotation genomic features}
#'
#'                  \item{\code{annotation}}{a named vector of percentages}
#'
#'                  \item{\code{precedence}}{a named vector of percentages}
#'
#'                  \item{\code{num.hierarchica}}{vector}
#'
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#' }
#' @name annotationByFeature-class
#' @rdname annotationByFeature-class
#' @export
setClass("annotationByFeature", 
			representation(members  = "matrix",
						   annotation = "numeric",
						   precedence = "numeric",
						   num.annotation = "numeric",
						   num.precedence = "numeric",
						   no.of.OlapFeat = "numeric",
						   perc.of.OlapFeat = "numeric"))

# ----------------------------------------------------------------------------------------------- #
#' An S4 class that information on overlap of target features with annotation features  
#'
#' This object is desgined to hold statistics and information about genomic feature overlaps
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target features with annotation genomic features}
#'
#'                  \item{\code{annotation}}{a named vector of percentages}
#'
#'                  \item{\code{precedence}}{a named vector of percentages}
#'
#'                  \item{\code{num.hierarchica}}{vector}
#'
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#'
#'                  \item{dist.to.TSS}{a data frame showing distances to TSS and gene/TSS names and strand}
#' }
#' @name annotationByGenicParts-class
#' @rdname annotationByGenicParts-class
#' @export
setClass("annotationByGenicParts", 
			representation(dist.to.TSS = "data.frame"), contains = "annotationByFeature")


#new.obj=new("annotationByGenicParts",
#            members=matrix(c(1,2,3,4)),annotation=c(1,2,0,3,4),
#            precedence=c(a=1,b=2,c=0,d=3,e=4),
#            num.annotation  =c(1,2,0,3,4),
#            num.precedence=c(1,2,0,3,4),
#            no.of.OlapFeat  =c(1,2,0,3,4),
#            perc.of.OlapFeat=c(1,2,0,3,4),
#            dist.to.TSS     =c(1,2,0,3,4) )

#' @rdname show-methods
#' @aliases show,annotationByGenicParts-method
setMethod("show", "annotationByGenicParts", 
			function(object) {

				cat("Summary of target set annotation with genic parts\n");
				cat("Rows in target set:", nrow(object@members), "\n")
				cat("-----------------------\n\n")
				
				cat("percentage of target features overlapping with annotation :\n")
				print(object@annotation)
				cat("\n\n")
				
				cat("percentage of target features overlapping with annotation\n")
				cat("(with promoter > exon > intron precedence) :\n"); 
				print(object@precedence)
				cat("\n\n")
				
				cat("percentage of annotation boundaries with feature overlap :\n")
				print(object@perc.of.OlapFeat);
				cat("\n\n")  
				
				cat("summary of distances to the nearest TSS :\n")
				print(summary(abs(object@dist.to.TSS[,2])))
				cat("\n")
})


#' @rdname show-methods
#' @aliases show,annotationByFeature-method
setMethod("show", "annotationByFeature", 
			function(object) {

			cat("summary of target set annotation with feature annotation\n")
			cat(nrow(object@members))
			cat(" rows in target set\n--------------\n")
			cat("--------------\n")
			
			cat("percentage of target features overlapping with annotation :\n")
			print(object@annotation)
			cat("\n\n")
			
			cat("percentage of annotation boundaries with feature overlap :\n")
			print(object@perc.of.OlapFeat)
			cat("\n\n")  
})


#######################################
# SECTION 2: S3 FUNCTIONS 
# these shouldn't be exported
#######################################

# ----------------------------------------------------------------------------------------------- #

annotatGrWithGeneParts <- function(gr, prom, exon, intron, strand=F){

	if( ! strand){strand(gr)="*"}
	memb = data.frame(matrix(rep(0,length(gr)*3),ncol=3) )
	colnames(memb)=c("prom","exon","intron")
	memb[countOverlaps(gr,prom) > 0,1] = 1
	memb[countOverlaps(gr,exon) > 0,2] = 1
	memb[countOverlaps(gr,intron) > 0,3] = 1

	# percentage of overlap for each region
	annotation = c( promoter 	= 100*sum(memb$prom>0)/nrow(memb),
					exon		= 100*sum(memb$exon>0)/nrow(memb),
					intron		= 100*sum(memb$intron>0)/nrow(memb),
					intergenic	= 100*sum(rowSums(memb)==0)/nrow(memb) )

	
	num.annotation = c( promoter	= sum(memb$prom>0),
						exon		= sum(memb$exon>0),
						intron		= sum(memb$intron>0),
						intergenic	= sum(rowSums(memb)==0) )

	precedence = c( promoter  = 100*sum(memb$prom>0)/nrow(memb),
					exon      = 100*sum(memb$exon>0 & memb$prom==0)/nrow(memb),
					intron    = 100*sum(memb$intron>0 & memb$exon==0 & memb$prom==0)/nrow(memb),
					intergenic= 100*sum(rowSums(memb)==0)/nrow(memb) )

	num.precedence = c( promoter  = sum(memb$prom>0),
						exon      = sum(memb$exon>0 & memb$prom==0),
						intron    = sum(memb$intron>0 & memb$exon==0 & memb$prom==0),
						intergenic= sum(rowSums(memb)==0) )  

	numberOfOlapFeat = c( promoter=sum(countOverlaps(prom,gr)>0),
						  exon=sum(countOverlaps(exon,gr)>0),
						  intron=sum(countOverlaps(intron,gr)>0) )
	percOfOlapFeat =100*numberOfOlapFeat/c(length(prom),length(exon),length(intron))
	
	return(list(members = memb,
				annotation=annotation,
				precedence=precedence,
				num.annotation=num.annotation,
				num.precedence=num.precedence,
				numberOfOlapFeat=numberOfOlapFeat,percOfOlapFeat=percOfOlapFeat))
}

# ----------------------------------------------------------------------------------------------- #
distance2NearestFeature<-function(g.idh,tss){

	elementMetadata(g.idh) = DataFrame(elementMetadata(g.idh),orig.row=1:length(g.idh))
	# get the row number column
	id.col = ncol( elementMetadata(g.idh))+5 

	# get the id column for tss in the merged data.frame below
	tss.id.col = ncol( elementMetadata(g.idh))+5+7 
	
	# get the id column for tss
	strnd.col=ncol( elementMetadata(g.idh))+5+7-2 

	# a merged data.frame is returned by this function
	met.tss = .nearest.2bed(g.idh, tss) 

	dist.col=ncol(met.tss)
	met.tss[met.tss$end < met.tss$start.y & met.tss$strand.y=="+", dist.col] = -1* met.tss[met.tss$end<met.tss$start.y & met.tss$strand.y=="+",dist.col]

	met.tss[met.tss$end > met.tss$start.y & met.tss$strand.y=="-",dist.col] = -1* met.tss[met.tss$end>met.tss$start.y & met.tss$strand.y=="-",dist.col]

	res = met.tss[order(met.tss[,id.col]),c(id.col,dist.col,tss.id.col,strnd.col)]
	names(res)=c("target.row", "dist.to.feature",  "feature.name", "feature.strand")

	return(res)
}


# ----------------------------------------------------------------------------------------------- #
# get the nearest features between a subject and query bed file
# get grange objects in bed zormat
# g.bed is GenomicRanges object, biatch!!!
# subject is also GenomicRanges object
# I think this function can be easily simplified
.nearest.2bed<-function(g.bed,subject){
	
	chrs1 = IRanges::levels(seqnames(g.bed))
	chrs2 = IRanges::levels(seqnames(subject))
	chrs  = intersect(chrs1, chrs2)
	res.df=NULL
	
	for(i in 1:length(chrs)){
		
		# find the nearest for this range
		ind = nearest(ranges(g.bed[seqnames(g.bed)==chrs[i],]),ranges(subject[seqnames(subject)==chrs[i],]))
		sbdf1 = IRanges::as.data.frame(g.bed[seqnames(g.bed)==chrs[i],] )
		sbdf2 = IRanges::as.data.frame(subject[seqnames(subject)==chrs[i],] )
		sbdf2 = sbdf2[ind,]
		names(sbdf2)=paste(names(sbdf2),".y",sep="")
		temp  = cbind(sbdf1,sbdf2)
		res.df = rbind(res.df,temp)
	}

	res.dist = rep(0,nrow(res.df))
	dist1 = abs(res.df$start-res.df$end.y)+1
	dist2 = abs(res.df$start.y-res.df$end)+1
	
	#total length
	totti = res.df$width + res.df$width.y 
	res.dist[pmax(dist1,dist2)>totti]=pmin(dist1,dist2)[pmax(dist1, dist2) > totti] 
	res.df = cbind(res.df,dist=res.dist)
	return(res.df)
}





#######################################
# SECTION 2: S4 FUNCTIONS 
#######################################


# ----------------------------------------------------------------------------------------------- #
#' Function to annotate given GRanges object with promoter, exon, intron & intergenic ranges
#'
#' @param target: A GRanges object storing chromosome locations to be annotated (e.g. chipseq peaks)
#' @param GRangesList.obj: A GRangesList object containing GRanges object for promoter, exons, introns and TSSes, or simply output of readTranscriptFeatures function
#' @param strand: If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FALSE)
#' @usage annotateWithGeneParts(target,GRangesList.obj,strand=F)
#' @return \code{annotationByGenicParts} object
#' 
#' @export
#' @docType methods
#' @rdname annotateWithGeneParts-methods
setGeneric("annotateWithGeneParts", function(target,GRangesList.obj,strand=F) standardGeneric("annotateWithGeneParts") )

#' @aliases annotateWithGeneParts,GRanges,GRangesList-method
#' @rdname annotateWithGeneParts-methods
setMethod("annotateWithGeneParts", 
		  signature(target = "GRanges", GRangesList.obj = "GRangesList"),
		  function(target, GRangesList.obj,strand){

			a.list = annotatGrWithGeneParts(target, 
												GRangesList.obj$promoters,
												GRangesList.obj$exons,
												GRangesList.obj$introns,
												strand=strand)
			dist2TSS = distance2NearestFeature(target,GRangesList.obj$TSSes)

			new("annotationByGenicParts",
				members 		= as.matrix(a.list$members),
				annotation      = a.list$annotation,
				precedence		= a.list$precedence,
				num.annotation  = a.list$num.annotation,
				num.precedence	= a.list$num.precedence,
				no.of.OlapFeat  = a.list$numberOfOlapFeat,
				perc.of.OlapFeat= a.list$percOfOlapFeat,
				dist.to.TSS     = dist2TSS )
})


# ----------------------------------------------------------------------------------------------- #
#'Given a GRangesList object it annotates each Range with gene annotation
#'
#' @param target: A GRangesList object storing chromosome locations to be annotated (e.g. chipseq peaks from multiple experiments)
#' @param GRangesList.obj: A GRangesList object containing GRanges object for promoter, exons, introns and TSSes, or simply output of readTranscriptFeatures function
#' @param strand: If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FALSE)
#'
#' @usage annotateWithGeneParts(target, GRangesList.obj, strand=F)
#' @return \code{annotationByGenicParts} object
#' @export
#'
#' @docType methods
#' @rdname annotateWithGeneParts-methods
setMethod("annotateWithGeneParts",
		  signature(target = "GRangesList", GRangesList.obj= "GRangesList"),
		  function(target, GRangesList.obj, strand){
		  
			l = lapply(target, function(x)annotateWithGeneParts(target, GRangesList.obj, strand))
			return(l)
})




# ----------------------------------------------------------------------------------------------- #
#' Function to annotate a given GRanges object with promoter,exon,intron & intergenic values
#'  
#' @param target    a granges object storing chromosome locations to be annotated
#' @param feature   a granges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param flank     a granges object storing chromosome locations of the flanks of the feature
#' @param feature.name     string for the name of the feature
#' @param flank.name     string for the name of the flanks
#' @param strand   If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @usage annotateWithFeatureFlank(target,feature,flank,feature.name="feat",flank.name="flank",strand=FALSE)
#' @return returns an \code{annotationByFeature} object
#' 
#' @export
#' @docType methods
#' @rdname annotateWithFeatureFlank-methods
setGeneric("annotateWithFeatureFlank", function(target,feature,flank,feature.name="feat",flank.name="flank",strand=FALSE) standardGeneric("annotateWithFeatureFlank") )

#' @aliases annotateWithFeatureFlank,GRanges,GRanges,GRanges-method
#' @rdname annotateWithFeatureFlank-methods
setMethod( "annotateWithFeatureFlank", 
			signature(target = "GRanges",feature="GRanges",flank="GRanges"),
			function(target, feature, flank,feature.name,flank.name,strand){

				if( ! strand )
					strand(target) = "*"
					memb=data.frame(matrix(rep(0,length(target)*2),ncol=2) )
					colnames(memb) = c(feature.name,flank.name)
					memb[countOverlaps(target,feature)>0,1]=1
					memb[countOverlaps(target,flank)>0,2]=1

					annotation=c(100*sum(memb[,1]>0)/nrow(memb),
								 100*sum(memb[,2]>0)/nrow(memb),
								 100*sum(rowSums(memb)==0)/nrow(memb) )
					names(annotation)=c(feature.name, flank.name, "other")

					num.annotation=c( sum(memb[,1]>0), 
									  sum(memb[,2]>0), 
									  sum(rowSums(memb)==0) )
					names(num.annotation)=c(feature.name,flank.name,"other")                      

					precedence=c(100*sum(memb[,1]>0)/nrow(memb) ,
								 100*sum(memb[,2]>0 & memb[,1]==0)/nrow(memb) ,
								 100*sum(rowSums(memb)==0)/nrow(memb) )
					names(precedence)=c(feature.name,flank.name,"other")

					num.precedence=c( sum(memb[,1]>0), 
									  sum(memb[,2]>0 & memb[,1]==0), 
									  sum(rowSums(memb)==0)  )
					names(num.precedence)=c(feature.name,flank.name,"other")                      

					numberOfOlapFeat = c(sum(countOverlaps(feature,target)>0),
					sum(countOverlaps(flank,target)>0) )
					names(numberOfOlapFeat) = c(feature.name,flank.name)
					percOfOlapFeat = 100*numberOfOlapFeat/c(length(feature),length(flank) )

					new("annotationByFeature",
						members         = as.matrix(memb),
						annotation      = annotation,
						precedence		= precedence,
						num.annotation  = num.annotation,
						num.precedence	= num.precedence,
						no.of.OlapFeat  = numberOfOlapFeat,
						perc.of.OlapFeat= percOfOlapFeat)

})


# ----------------------------------------------------------------------------------------------- #
#' Function to annotate given GRanges object with a given genomic feature
#' 
#' @param target   a GRanges object storing chromosome locations to be annotated
#' @param feature  a GRanges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param strand   If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @param extend   specifiying a positive value will extend the feature on both sides as much as \code{extend}
#' @param feature.name name of the annotation feature. For example: H3K4me1,CpGisland etc.
#' @usage annotateWithFeature(target,feature,strand=FALSE,extend=0,feature.name="feat1")
#' @return returns an \code{annotationByFeature} object
#' 
#' @export
#' @docType methods
#' @rdname annotateWithFeature-methods
setGeneric("annotateWithFeature", function(target,feature,strand=FALSE,extend=0,feature.name="feat1") standardGeneric("annotateWithFeature") )

#' @aliases annotateWithFeature,GRanges,GRanges-method
#' @rdname annotateWithFeature-methods
setMethod("annotateWithFeature", 
		   signature(target = "GRanges",feature="GRanges"),
		   function(target, feature, strand,extend,feature.name){

				if( ! strand){strand(target)="*"}
					memb=rep(0,length(target))

				if(extend>0){
					start(feature) = start(feature)- extend
					end(feature)   = end(feature)  + extend
				}
				memb[countOverlaps(target,feature) > 0] = 1

				annotation = c( 100*sum(memb >  0)/length(memb) ,
								100*sum(memb == 0)/length(memb) )
				
				num.annotation = c( sum( memb >  0),
									sum( memb == 0) )
				names(annotation) = c(feature.name, "other")

				numberOfOlapFeat = c(sum(countOverlaps(feature,target)>0))
				percOfOlapFeat = 100*numberOfOlapFeat/c(length(feature))

				new("annotationByFeature",
					members         = as.matrix(memb),
					annotation      = annotation,
					precedence		= annotation,
					num.annotation  = num.annotation,
					num.precedence	= num.annotation,
					no.of.OlapFeat  = numberOfOlapFeat,
					perc.of.OlapFeat= percOfOlapFeat)
})

# ACCESSOR FUNCTIONS
#annotationByFeature
#annotationBygenicparts


# ----------------------------------------------------------------------------------------------- #
#' Get the membership slot of annotationByFeature
#'
#' Membership slot defines the overlap of target features with annotation features
#'  For example, if a target feature overlaps with an exon
#' 
#' @param x a \code{annotationByFeature}  object
#' 
#' @return RETURNS a matrix showing overlap of target features with annotation features. 1 for overlap, 0 for non-overlap
#' 
#' @usage getMembers(x)
#'
#' @export
#' @docType methods
#' @rdname getMembers-methods                      
setGeneric("getMembers", def=function(x) standardGeneric("getMembers"))

#' @aliases getMembers,annotationByFeature-method
#' @rdname getMembers-methods
setMethod("getMembers", 
			signature(x = "annotationByFeature"),
			function(x){
				return(x@members)
})

# ----------------------------------------------------------------------------------------------- #
#' Get the percentage of target features overlapping with annotation from annotationByFeature
#'
#' This function retrieves percentage/number of target features overlapping with annotation
#'  
#' @param x a \code{annotationByFeature}  object
#' @param percentage TRUE|FALSE. If TRUE percentage of target features will be returned. If FALSE, number of target features will be returned
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation features when calculating numbers (with promoter>exon>intron precedence)
#' That means if a feature overlaps with a promoter it will be counted as promoter overlapping only, or if it is overlapping with a an exon but not a promoter, 
#' it will be counted as exon overlapping only whether or not it overlaps with an intron.
#'
#' @usage getTargetAnnotationStats(x,percentage=TRUE,precedence=TRUE)
#'
#' @return RETURNS  a vector of percentages or counts showing quantity of target features overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname getTargetAnnotationStats-methods
setGeneric("getTargetAnnotationStats", def=function(x,percentage=TRUE,precedence=TRUE) standardGeneric("getTargetAnnotationStats"))

#' @rdname getTargetAnnotationStats-methods
#' @aliases getTargetAnnotationStats,annotationByFeature-method
setMethod("getTargetAnnotationStats", 
			signature(x = "annotationByFeature"),
			function(x,percentage ,precedence ){                      
			
			if(percentage){
				if(precedence){
					return(x@precedence)
				}else{
					return(x@annotation)
				}
			}else{
				if(precedence){
					return(x@num.precedence)
				}else{
					return(x@num.annotation)
				}
			}
})


# ----------------------------------------------------------------------------------------------- #
#' Get the percentage/count of annotation features overlapping with target features from annotationByFeature
#'
#' This function retrieves percentage/number of annotation features overlapping with targets. 
#' For example, if \code{annotationByFeature}  object is containing statistics of differentially methylated 
#' regions overlapping with gene annotation. This function will return number/percentage of introns,exons and promoters
#' overlapping with differentially methylated regions.
#'  
#' @param x a \code{annotationByFeature}  object
#' @param percentage TRUE|FALSE. If TRUE percentage of annotation features will be returned. If FALSE, number of annotation features will be returned
#'
#' @return RETURNS  a vector of percentages or counts showing quantity of annotation features overlapping with target features
#' @usage getFeatsWithTargetsStats(x,percentage=TRUE)
#' @export
#' @docType methods
#' @rdname getFeatsWithTargetsStats-methods
setGeneric("getFeatsWithTargetsStats", def=function(x,percentage=TRUE) standardGeneric("getFeatsWithTargetsStats"))


#' @rdname getFeatsWithTargetsStats-methods
#' @aliases getFeatsWithTargetsStats,annotationByFeature-method
setMethod("getFeatsWithTargetsStats", 
			signature(x = "annotationByFeature" ),
			function( x,percentage ){                      
			
				if(percentage){
					return(x@perc.of.OlapFeat)
				}else{
					return(x@no.of.OlapFeat)                        
				}
})


# ----------------------------------------------------------------------------------------------- #
#' Get distance to nearest TSS and gene id from annotationByGenicParts
#'
#' This accessor function gets the nearest TSS, its distance to target feature,strand and name of TSS/gene from annotationByGenicParts object
#' @param x a \code{annotationByGenicParts}  object
#' 
#' @return RETURNS a data.frame containing row number of the target features,distance of target to nearest TSS, TSS/Gene name, TSS strand
#' @usage getAssociationWithTSS(x)
#' @aliases getAssociationWithTSS,-methods getAssociationWithTSS,annotationByGenicParts-method
#' @export
#' @docType methods
#' @rdname annotationByGenicParts-methods
setGeneric("getAssociationWithTSS", def=function(x) standardGeneric("getAssociationWithTSS"))

#' @rdname annotationByGenicParts-methods
#' @docType methods
#' @aliases getAssociationWithTSS,annotationByGenicParts-method
setMethod("getAssociationWithTSS", 
			signature(x = "annotationByGenicParts"),
			function(x){
				return(x@dist.to.TSS)
})



# PLOTTING FUNCTIONS

# ----------------------------------------------------------------------------------------------- #
#' Plot annotation categories from annotationByGenicParts or annotationByFeature
#'
#' This function plots a pie or bar chart for showing percentages of targets annotated by genic parts or other query features
#' @param x a \code{annotationByFeature} or  \code{annotationByGenicParts} object
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation features when calculating numbers (with promoter>exon>intron precedence). 
#'  This option is only valid when x is a \code{annotationByGenicParts} object
#' @param col a vector of colors for piechart or the bar plot
#' @param ... graphical parameters to be passed to \code{pie} or \code{barplot} functions
#'
#' usage  plotTargetAnnotation(x,precedence=TRUE,col,...)
#'
#' @return plots a piechart or a barplot for percentage of the target features overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname plotTargetAnnotation-methods
setGeneric("plotTargetAnnotation", def=function(x,precedence=TRUE,col=rainbow(length(x@annotation)),...) standardGeneric("plotTargetAnnotation"))

#' @rdname plotTargetAnnotation-methods
#' @docType methods
#' @aliases plotTargetAnnotation,annotationByFeature-method
setMethod("plotTargetAnnotation", 
			signature(x = "annotationByFeature"),
			function(x,precedence,col,...){
				
				props=getTargetAnnotationStats(x,precedence)

			if(precedence | ( is(x,"annotationByFeature") & !is(x,"annotationByGenicParts")) ){
				slice.names=names(props)
				names(props)=paste( paste(round(props),"%"),sep=" ")

				pie(props,cex=0.9,col=col,...)
				legend("topright",legend=slice.names,fill=col )

			}else{

				mids = barplot(props,col=col,...)  
				text(mids, y = props,labels=paste(paste(round(props),"%",sep="")),pos=1) 
			}
})


# ----------------------------------------------------------------------------------------------- #
#' Plots the enrichment of each feature in the set in the gene annotation
#'
#' This function plots a heatmap of enrichment of each range in given gene feature
#' @param l a \code{list} of  \code{annotationByGenicParts} objects
#' @param cluster TRUE/FALSE. If TRUE the heatmap is going to be clustered using hierarchical clustering
#' @param col a vector of two colors that will be used for interpolation. The first color is the lowest one, the second is the highest one

#' usage  plotGenicAnnotation(l, cluster=FALSE, col=c('white','cornflowerblue'),...)
#' 
#' @export
#' @docType methods
#' @rdname plotGenicAnnotation-methods
setGeneric("plotGenicAnnotation", def=function(l, cluster=FALSE, col=c('white','cornflowerblue'), ...)standardGeneric("plotGenicAnnotation"))

#' @rdname plotGenicAnnotation-methods
#' @docType methods
#' @aliases plotGenicAnnotation,annotationByFeature-method
setMethod("plotGenicAnnotation", 
			signature(l="list"),
			function(l, cluster, col, ...){
			
				if(!all(unlist(lapply(l, class)) == "annotationByGenicParts"))
					stop("All elements of the input list need to be annotationByGenicParts-class")
			
				d = do.call(rbind, lapply(l, function(x)x@precedence))
				
				ind = 1:nrow(d)
				if(cluster == TRUE){
					h = hclust(dist(d))
					ind = h$order
				}
				m = reshape2::melt(data.frame(d))
				m = melt(d)
				
				p <- ggplot(m, aes(x=X2, y=X1, fill=value, colour="white")) + 
					 scale_fill_gradient(low =col[1], high = col[2]) + 
					 scale_y_discrete(limits = rownames(d)[ind] ) + 
					 opts(axis.title.x=theme_text(colour='white'), 
						  axis.text.x=theme_text(colour='black', face='bold'), 
						  axis.text.y=theme_text(colour='black'), 
						  axis.title.y=theme_text(colour='white', face='bold'))
				p + geom_tile(color='white') 
})

# SECTION 3:
# annotate ML objects with annotations read-in and converted to GRanges objects


