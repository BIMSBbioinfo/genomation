# PART THAT DEALS with annotation of genomic features

##############################################################################
# SECTION 1:
# reading annotation to GRanges
# makes GRanges object from a given bed6 or bed12 file to granges object
##############################################################################

#######################################
# SECTION 1: S3 functions
#######################################

# ---------------------------------------------------------------------------- #
# extracts exons from a bed12 file and puts them into GRanges object
bed12ToExons<-function(ref){

	ref=unique(ref)

	# apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
	b.start.size = cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),
                       as.integer(unlist(strsplit(as.character(ref$V11),"," ))))
	rep.ref = ref[rep(1:nrow(ref),ref[,10]),] # replicate rows occurs as many as its exon number


	exon.id = unlist( mapply( function(x,y)
                              if(x=="+"){return(1:y)}
                              else{return(y:1)}, 
                              ref[,6],ref[,10]))
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

# ---------------------------------------------------------------------------- #
# extracts introns from a bed12 file and puts them into GRanges object
bed12ToIntrons<-function(ref){

	#remove the genes with one exon only (they won't have any introns)
	ref=ref[ref[,10]>1,]
	ids=paste(ref[,1],ref[,2],ref[,3],ref[,4],sep="")
	ref=cbind(ref,id=ids)
	ref=unique(ref)

	# apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
	b.start.size=cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),
                     as.integer(unlist(strsplit(as.character(ref$V11),"," ))))

	# replicate rows occurs as many as its exon number
	rep.ref = ref[rep(1:nrow(ref),ref[,10]),] 

	exon.id = unlist( mapply( function(x,y) 
                                if(x=="+"){return(1:y)}
                                else{return(y:1)} ,
                                ref[,6],ref[,10]  ) )
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

# ---------------------------------------------------------------------------- #
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


# ---------------------------------------------------------------------------- #
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
				grange = GRanges(seqnames=bed$V1,
                         ranges=IRanges(start=bed$V2+1, end=bed$V3),
                         strand=rep("*",nrow(bed)), score=bed$V5,name=bed$V4 )
			}
			if(ncol(bed)==4){
				grange = GRanges(seqnames=bed$V1,
                         ranges=IRanges(start=bed$V2+1, end=bed$V3),
                         strand=rep("*",nrow(bed)),name=bed$V4)
			}    
			if(ncol(bed)==3){
				grange = GRanges(seqnames=bed$V1,
                         ranges=IRanges(start=bed$V2+1, end=bed$V3),
                         strand=rep("*",nrow(bed)))
			}                           
			return(grange)
})

# ---------------------------------------------------------------------------- #
#' convert a data frame read-in from a bed file to a GRanges object for exons
#'  
#' @param bed.df a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convertBed2Exons(bed.df)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname convertBed2Exons-methods
setGeneric("convertBed2Exons",
                function(bed.df) standardGeneric("convertBed2Exons"))

#' @aliases convertBed2Exons,data.frame-method
#' @rdname convertBed2Exons-methods
setMethod("convertBed2Exons" ,
			signature(bed.df = "data.frame" ),
			function(bed.df){

			if(! checkBedValidity(bed.df,"exon"))
				stop("this is not a valid bed file")
				
			bed12ToExons(bed.df)
})


# ---------------------------------------------------------------------------- #
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
setGeneric("convertBed2Introns",
                  function(bed.df) standardGeneric("convertBed2Introns"))

#' @aliases convertBed2Introns,data.frame-method
#' @rdname convertBed2Introns-methods
setMethod("convertBed2Introns",
			signature(bed.df = "data.frame" ),
			function(bed.df){

			if(! checkBedValidity(bed.df,"exon"))
				stop("this is not a valid bed file")
				
			bed12ToIntrons(bed.df)
})



# ---------------------------------------------------------------------------- #
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
setGeneric("getFlanks", 
                function(grange,flank=2000,clean=T) 
                          standardGeneric("getFlanks"))

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




##############################################################################
# SECTION 2:
# annotate granges objects with annotations that read-in and converted to GRanges objects
##############################################################################

#' @rdname show-methods
#' @aliases show,annotationByGenicParts-method
setMethod("show", "annotationByGenicParts", 
			function(object) {

				message("Summary of target set annotation with genic parts\n");
				message("Rows in target set:", nrow(object@members), "\n")
				message("-----------------------\n\n")
				
				message("percentage of target features overlapping with annotation :\n")
				print(round(object@annotation,2))
				message("\n\n")
				
				message("percentage of target features overlapping with annotation\n")
				message("(with promoter > exon > intron precedence) :\n"); 
				print(round(object@precedence,2))
				message("\n\n")
				
				message("percentage of annotation boundaries with feature overlap :\n")
				print(round(object@perc.of.OlapFeat,2))
				message("\n\n")  
				
				message("summary of distances to the nearest TSS :\n")
				print(summary(abs(object@dist.to.TSS[,2])))
				message("\n")
})


#' @rdname show-methods
#' @aliases show,annotationByFeature-method
setMethod("show", "annotationByFeature", 
			function(object) {

			message("summary of target set annotation with feature annotation\n")
			message(nrow(object@members))
			message(" rows in target set\n--------------\n")
			message("--------------\n")
			
			message("percentage of target features overlapping with annotation :\n")
			print(round(object@annotation,2))
			message("\n\n")
			
			message("percentage of annotation boundaries with feature overlap :\n")
			print(round(object@perc.of.OlapFeat))
			message("\n\n")  
})


#######################################
# SECTION 2: S3 FUNCTIONS 
# these shouldn't be exported
#######################################

# ---------------------------------------------------------------------------- #

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

# ---------------------------------------------------------------------------- #
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
	met.tss[met.tss$end < met.tss$start.y & met.tss$strand.y=="+", dist.col] = 
    -1* met.tss[met.tss$end<met.tss$start.y & met.tss$strand.y=="+",dist.col]

	met.tss[met.tss$end > met.tss$start.y & met.tss$strand.y=="-",dist.col] = 
    -1* met.tss[met.tss$end>met.tss$start.y & met.tss$strand.y=="-",dist.col]

	res = met.tss[order(met.tss[,id.col]),c(id.col,dist.col,tss.id.col,strnd.col)]
	names(res)=c("target.row", "dist.to.feature",  "feature.name", "feature.strand")

	return(res)
}


# ---------------------------------------------------------------------------- #
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
		ind = nearest(ranges(g.bed[seqnames(g.bed)==chrs[i],]),
                  ranges(subject[seqnames(subject)==chrs[i],]))
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


# ---------------------------------------------------------------------------- #
#' Annotate given  object with promoter, exon, intron and intergenic regions
#' 
#' The function annotates \code{GRangesList} or \code{GRanges} object as 
#' overlapping with promoter,exon,intron or intergenic regions.
#'
#' @param target  \code{\link{GRanges}} or \code{\link{GRangesList}}
#'                object storing chromosome locations
#'                to be annotated (e.g. chipseq peaks)
#' @param feature GRangesList object containing GRanges object for
#'                promoter, exons, introns and transcription start sites, 
#'                or simply output of readTranscriptFeatures function
#' @param strand If set to TRUE, annotation features and target features will be 
#'        overlapped based on strand  (def:FALSE)
#' @param intersect.chr boolean, whether to select only chromosomes that are 
#'        common to feature and target. FALSE by default
#' 
#' @return \code{annotationByGenicParts} object or a list of 
#'         \code{annotationByGenicParts}
#'         objects if target is a  \code{\link{GRangesList}} object.
#' 
#' @export
#' @docType methods
#' @rdname annotateWithGeneParts-methods
setGeneric("annotateWithGeneParts", 
                function(target,feature,strand=F, intersect.chr=FALSE)
                    standardGeneric("annotateWithGeneParts") )

#' @aliases annotateWithGeneParts,GRanges,GRangesList-method
#' @rdname annotateWithGeneParts-methods
setMethod("annotateWithGeneParts", 
		  signature(target = "GRanges", feature = "GRangesList"),
		  function(target, feature, strand, intersect.chr=FALSE){

		  if(intersect.chr){
		    message('intersecting chromosomes...')
		    chrs = intersect(seqnames(target), 
                         unique(unlist(lapply(feature), seqnames)))
        
		    if(length(chr) == 0)
		      stop('target and feature do not have intersecting chromosomes')
		    target=target[seqnames(target) %in% chrs]
		    feature = lapply(feature, function(x)x[seqnames(x) %in% chrs])
		  }  
        
			a.list = annotatGrWithGeneParts(target, 
			                                feature$promoters,
			                                feature$exons,
			                                feature$introns,
												              strand=strand)
			dist2TSS = distance2NearestFeature(target,feature$TSSes)

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



#' @aliases annotateWithGeneParts,GRangesList,GRangesList-method
#' @rdname annotateWithGeneParts-methods
setMethod("annotateWithGeneParts",
		  signature(target = "GRangesList", GRangesList.obj= "GRangesList"),
		  function(target, GRangesList.obj, strand){
		  
			l = lapply(target, function(x)
                          annotateWithGeneParts(target, GRangesList.obj, strand))
			return(l)
})




# ---------------------------------------------------------------------------- #
#' Function to annotate a given GRanges object with promoter,exon,intron & intergenic values
#'  
#' @param target    a granges object storing chromosome locations to be annotated
#' @param feature   a granges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param flank     a granges object storing chromosome locations of the flanks of the feature
#' @param feature.name     string for the name of the feature
#' @param flank.name     string for the name of the flanks
#' @param strand   If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @param intersect.chr boolean, whether to select only chromosomes that are common to feature and target. FALSE by default
#' @usage annotateWithFeatureFlank(target,feature,flank,feature.name="feat",flank.name="flank",strand=FALSE)
#' @return returns an \code{annotationByFeature} object
#' 
#' @export
#' @docType methods
#' @rdname annotateWithFeatureFlank-methods
setGeneric("annotateWithFeatureFlank", 
                          function(target,
                                   feature,
                                   flank,
                                   feature.name="feat",
                                   flank.name="flank",
                                   strand=FALSE) 
                            standardGeneric("annotateWithFeatureFlank") )

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


# ---------------------------------------------------------------------------- #
#' Function to annotate given GRanges object with a given genomic feature
#' 
#' @param target   a GRanges object storing chromosome locations to be annotated
#' @param feature  a GRanges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param strand   If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @param extend   specifiying a positive value will extend the feature on both sides as much as \code{extend}
#' @param feature.name name of the annotation feature. For example: H3K4me1,CpGisland etc.
#' @param intersect.chr boolean, whether to select only chromosomes that are common to feature and target. FALSE by default

#' @usage annotateWithFeature(target,feature,strand=FALSE,extend=0,feature.name="feat1")
#' @return returns an \code{annotationByFeature} object
#' 
#' 
#' @examples
#' data(cpgi)
#' data(promoters)
#' annot = annotateWithFeature(cpgi, promoters)
#' 
#' 
#' @export
#' @docType methods
#' @rdname annotateWithFeature-methods
setGeneric("annotateWithFeature", 
                    function(target,
                             feature,
                             strand=FALSE,
                             extend=0,
                             feature.name="feat1", 
                             intersect.chr=FALSE) 
                      standardGeneric("annotateWithFeature") )

#' @aliases annotateWithFeature,GRanges,GRanges-method
#' @rdname annotateWithFeature-methods
setMethod("annotateWithFeature", 
		   signature(target = "GRanges",feature="GRanges"),
		   function(target, feature, strand,extend,feature.name,intersect.chr){

				if( ! strand){strand(target)="*"}
					memb=rep(0,length(target))

				if(extend>0){
				  message('extending features...')
          feature = resize(feature, width=(width(feature)+2*extend), fix='center')
				}
        
        # selects common chromosomes for target and feature
        if(intersect.chr){
          message('intersecting chromosomes...')
          chrs = intersect(seqnames(feature), seqnames(target))
          if(length(chr) == 0)
            stop('target and feature do not have intersecting chromosomes')
          target=target[seqnames(target) %in% chrs]
          feature = feature[seqnames(feature) %in% chrs]
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


# ---------------------------------------------------------------------------- #
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

# ---------------------------------------------------------------------------- #
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
setGeneric("getTargetAnnotationStats", 
                  function(x,percentage=TRUE,precedence=TRUE) 
                    standardGeneric("getTargetAnnotationStats"))

#' @rdname getTargetAnnotationStats-methods
#' @aliases getTargetAnnotationStats,annotationByFeature-method
setMethod("getTargetAnnotationStats", 
			signature(x = "annotationByFeature"),
			function(x, percentage ,precedence ){                      
			
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


# ---------------------------------------------------------------------------- #
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
setGeneric("getFeatsWithTargetsStats", 
                  function(x,percentage=TRUE) 
                    standardGeneric("getFeatsWithTargetsStats"))


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


# ---------------------------------------------------------------------------- #
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
setGeneric("getAssociationWithTSS", 
                  function(x) 
                    standardGeneric("getAssociationWithTSS"))

#' @rdname annotationByGenicParts-methods
#' @docType methods
#' @aliases getAssociationWithTSS,annotationByGenicParts-method
setMethod("getAssociationWithTSS", 
			signature(x = "annotationByGenicParts"),
			function(x){
				return(x@dist.to.TSS)
})



# PLOTTING FUNCTIONS

# ---------------------------------------------------------------------------- #
#' Plot annotation categories from annotationByGenicParts or annotationByFeature
#'
#' This function plots a pie or bar chart for showing percentages of targets 
#' annotated by genic parts or other query features
#' 
#' @param x a \code{annotationByFeature} or  
#'        \code{annotationByGenicParts} object
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation 
#'        features when calculating numbers 
#'        (with promoter>exon>intron precedence). 
#'        This option is only valid when x is a 
#'        \code{annotationByGenicParts} object
#' @param col a vector of colors for piechart or the bar plot
#' @param ... graphical parameters to be passed to \code{pie} 
#'            or \code{barplot} functions
#'
#' usage  plotTargetAnnotation(x,precedence=TRUE,col,...)
#'
#' @return plots a piechart or a barplot for percentage of 
#'         the target features overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname plotTargetAnnotation-methods
setGeneric("plotTargetAnnotation", 
              function(x,
                       precedence=TRUE,
                       col=rainbow(length(x@annotation)),...) 
                       standardGeneric("plotTargetAnnotation"))

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


# ---------------------------------------------------------------------------- #
#' Plots the enrichment of each feature in the set in the gene annotation
#'
#' This function plots a heatmap of enrichment of each range in given gene feature
#' 
#' @param l a \code{list} of  \code{annotationByGenicParts} objects
#' @param cluster TRUE/FALSE. If TRUE the heatmap is going to be clustered 
#'        using hierarchical clustering
#' @param col a vector of two colors that will be used for interpolation. 
#'        The first color is the lowest one, the second is the highest one

#' usage  plotGeneAnnotation(l, cluster=FALSE, 
#'        col=c('white','cornflowerblue'),...)
#' 
#' @export
#' @docType methods
#' @rdname plotGeneAnnotation-methods
setGeneric("plotGeneAnnotation", 
                  function(l, 
                           cluster=FALSE, 
                           col=c('white','cornflowerblue'))
                    standardGeneric("plotGeneAnnotation"))

#' @rdname plotGeneAnnotation-methods
#' @docType methods
#' @aliases plotGeneAnnotation,annotationByFeature-method
setMethod("plotGeneAnnotation", 
			signature(l="list"),
			function(l, cluster, col){
			
				if(!all(unlist(lapply(l, class)) == "annotationByGenicParts"))
					stop("All elements of the input list need to be annotationByGenicParts-class")
			
				d = data.frame(do.call(rbind, lapply(l, function(x)x@precedence)))
				d = data.frame(SampleName=rownames(d), d)
        
				ind = 1:nrow(d)
				if(cluster == TRUE){
					h = hclust(dist(d[,-1]))
					ind = h$order
				}
        d = data.frame(ind, d)
				m = reshape2::melt(d, id.vars=c('ind', 'SampleName'))
        colnames(m)[3] = 'Position'
        
				p = ggplot(m, aes(x=Position, y=SampleName, fill=value, colour="white")) +
				  scale_fill_gradient(low =col[1], high = col[2]) +
					 theme(
              axis.title.x = element_text(colour='white'), 
						  axis.text.x  = element_text(colour='black', face='bold'), 
						  axis.text.y  = element_text(colour='black'), 
						  axis.title.y = element_text(colour='white', face='bold'))
				p = p + geom_tile(color='white')  
        print(p)
})

