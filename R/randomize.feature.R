# randomize the bed coordinates and calculate associations with a target

##############################################################################
## S3 functions to be used in S4 stuff
##############################################################################



##############################################################################
## S4
##############################################################################
#' function that randomizes the genomic coordinates
#'
#' This function randomly distributes the coordinates of genomic features which is stored in a \code{GRanges} object. The randomization can be constrained by supplied arguments.
#'
#' @param feature a GRanges object to be randomized
#' @param chrom.sizes sizes of chromosomes as a named vector (names are chromsomes names and elements of the vectors are lengths). If the vector of chromosome sizes isnot given seqlengts from the GRanges object will be used. If the seqlengths is not assigned, as the end of each chr the end last feature on each chr will be taken.
#' @param stranded if FALSE, all of the returned features will be strandless (will have "*" in the strand slot)
#' @param keep.strand.prop If TRUE strands will have the same proportion as the features
#' @param keep.chrom If TRUE, number of features and randomized features for a chromosome will match. Currently seeting this to FALSE is not supported.
#' @param exclude A GRanges object where no randomized feature should overlap, can be gaps or unmappable regions in the genome as an example.
#' @param include A GRanges object which defines the boundaries of randomized features
#' @param seed random number generator seed
#' @usage randomize.feature(feature, chrom.sizes=NULL, stranded=TRUE,  keep.strand.prop=TRUE, keep.chrom=TRUE, exclude=NULL, include=NULL, seed=NULL)
#' @return returns a GRanges object which is randomized version of the feature
#' @export
#' @docType methods
#' @rdname randomize.feature-methods
setGeneric("randomize.feature", function(feature,chrom.sizes=NULL,stranded=TRUE,keep.strand.prop=TRUE,keep.chrom=TRUE,exclude=NULL,include=NULL,seed=NULL) standardGeneric("randomize.feature") )

#' @aliases randomize.feature,GRanges-method
#' @rdname randomize.feature-methods
setMethod("randomize.feature", 
			signature(feature = "GRanges"),
			function( feature,chrom.sizes, stranded, keep.strand.prop, keep.chrom,exclude, include,seed ){

			if(!is.null(seed) & is.numeric(seed))
				set.seed(seed)

			# if there are chromosome sizes defined
			if(!is.null(chrom.sizes) ){

				max.inc = GRanges(seqnames = names(chrom.sizes),
								  ranges = IRanges( start=rep(1,length(chrom.sizes)), 
								  width=chrom.sizes))
			}


			# find the ranges of inclusion if include and/or exclude is given
			if( (  !is.null(include) & is(include,"GRanges") )  | ( !is.null(exclude) & is(exclude,"GRanges")  )){

				# define maximal inclusion from feature ranges or chromose.sizes
				if(!exists("max.inc")){
					chrs = as.character(unique(seqnames(feature)))
					my.starts=c()
					my.ends=c()
					for( i in 1:length(chrs) ){
						my.set	 = feature[seqnames(feature)==chrs[i], ]
						my.ends  = c(my.ends, max( end(my.set) ) )
						my.starts= c(my.starts, min( start(my.set) ) )
					}

					max.inc = GRanges(seqnames=chrs,
									  ranges=IRanges( start=rep(1,length(my.ends)), end=my.ends))
				}
				# if include is defined subtract intersect maximal inclusion with include ( or define include as maximal inclusion)
				if((  !is.null(include) & is(include,"GRanges") )  ){
					max.inc=intersect(max.inc,include)
				}

				# if exlcude is defined subtract this from maximal inclusion
				if((  !is.null(exclude) & is(exclude,"GRanges") )  ){
					max.inc=setdiff(max.inc,exclude)
				}
			}

			# if a maximal inclusion range NOT defined, find the range accordingly
			if(keep.chrom){
				if( ! exists("max.inc")){

					chrs=as.character(unique(seqnames(feature)))
					Ns=numeric(length(chrs));names(Ns)=chrs
					ends=numeric(length(chrs))

					# get the arguments to runif and construct GRanges after exiting the loop, costs less than appending
					for( i in 1:length(chrs) ){
						my.set= feature[seqnames(feature)==chrs[i], ]
						end.on.chr=max( end(my.set) )-max(width(my.set))
						Ns[i]=length(my.set)
						ends[i]=end.on.chr
					}
					
					my.runif=function(n,y) round(runif(n, min=1, max=y))

					# if you want to keep the strand proportion as well
					
					if(keep.strand.prop &  stranded){
						randomized = GRanges(seqnames = rep(names(Ns),Ns),
											 ranges = IRanges( start=unlist(mapply(my.runif,Ns,ends)), width  = width(feature)),
											 strand = strand(feature),
											 elementMetadata(feature))
					
					# if you  want to randomly assign strands
					}else if(!stranded){ 

						randomized=GRanges(seqnames	= rep(names(Ns),Ns),
										   ranges	= IRanges( start=unlist(mapply(my.runif,Ns,ends)), width	= width(feature)),
										   strand	= "*",
										   elementMetadata(feature))

					}else{
						
						randomized=GRanges(seqnames= rep(names(Ns),Ns),
										   ranges  = IRanges( start=unlist(mapply(my.runif,Ns,ends)), width   = width(feature)),
										   strand  = sample(unique(strand(feature)),sum(Ns),replace=T),
										   elementMetadata(feature))
					}
					return(randomized)
				
				# if a maximal inclusion range defined, find the range accordingly
				}else{ 

					 # get only chrs that are in both constraint and feature
					chrs1=as.character(unique(seqnames(max.inc)))
					chrs2=as.character(unique(seqnames(feature)))
					chrs=intersect(chrs1, chrs2)

					max.inc1 = max.inc[seqnames(max.inc) %in% chrs,]
					feature1 = feature[seqnames(feature) %in% chrs,]

					# create that will be feed to runif for random number generation
					Ns=numeric(length(max.inc1));names(Ns)=seqnames(max.inc1)
					ends=numeric(length(max.inc1))
					starts=numeric(length(max.inc1))
					my.index=1
					for( i in 1:length(chrs) ){

						Ns.chr= sum(seqnames(feature)==chrs[i]) # number of features on chr
						my.widths= width(max.inc[seqnames(max.inc)==chrs[i], ]) # widths of constraints on chr
						Ns.sub=round(Ns.chr*(my.widths/sum(my.widths))) # number of adjusted features per constraint range
						my.starts=start(max.inc[seqnames(max.inc)==chrs[i], ]) # starts of constraint on chr
						my.ends=end(max.inc[seqnames(max.inc)==chrs[i], ]) # end of constraints on chr

						if(sum(Ns.sub) > (Ns.chr))
							Ns.sub[which.max(Ns.sub)]=Ns.sub[which.max(Ns.sub)]-(sum(Ns.sub)-(Ns.chr))

						if(sum(Ns.sub) < (Ns.chr))
							Ns.sub[which.max(Ns.sub)]=Ns.sub[which.max(Ns.sub)]+( (Ns.chr)-sum(Ns.sub) )

						# update vectors that will be fed to runif
						Ns[my.index:(my.index+length(my.widths)-1)]=Ns.sub
						starts[my.index:(my.index+length(my.widths)-1)]=my.starts
						ends[my.index:(my.index+length(my.widths)-1)]=my.ends
						my.index=my.index+length(my.widths)
					}

					# remove ranges with no counts
					# since we distribute number of elements on ranges proportional to their length
					# short ranges may not have any elements on them
					starts=starts[Ns >0]
					ends  =ends[Ns>0]
					Ns    =Ns[Ns>0]

					my.runif=function(n,x,y) round(runif(n, min=x, max=y))

					# if you want to keep the strand proportion as well
					#
					if(keep.strand.prop & stranded){
						randomized = GRanges(seqnames=rep(names(Ns),Ns),
										   ranges  =IRanges( start=unlist(mapply(my.runif,Ns,starts,ends)), 
										   width= width(feature1)),
										   strand  =strand(feature1),
										   elementMetadata(feature1))
					 
					 # if you  want to randomly assign strands
					}else if(! stranded){
						randomized = GRanges(seqnames=rep(names(Ns),Ns),
											 ranges  = IRanges( start=unlist(mapply(my.runif,Ns,starts,ends)), 
											 width= width(feature1)),
											 strand="*",
											 elementMetadata(feature1))
					}else{

						randomized = GRanges(seqnames = rep(names(Ns),Ns),
											 ranges   = IRanges( start=unlist(mapply(my.runif,Ns,starts,ends)), 
											 width = width(feature1)),
											 strand  = sample(unique(strand(feature1)),sum(Ns),replace=T),
											 elementMetadata(feature1))
					}
					return(randomized)
				}

			}else{
				stop("UNSUPORRTED OPTION: 'keep.chrom' option has to be TRUE, FALSE is not supported yet.")
			}

})

