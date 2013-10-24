
# ------------------------------------------------------------------------------------------------ #
# randomizeFeature
#' function that randomizes the genomic coordinates
#'
#' This function randomly distributes the coordinates of genomic features which
#' is stored in a \code{GRanges} object. The randomization can be constrained by
#'  supplied arguments.
#' The function is still in Beta mode - the regions can overlap excluded regions, and the randomized regions are not disjoint. Please take care that the excluded and included regions are not too strict when compared to the total width of the ranges.
#'
#'
#' @param feature a GRanges object to be randomized
#' @param chrom.sizes sizes of chromosomes as a named vector (names are
#'        chromsomes names and elements of the vectors are lengths). , if not
#'        given sizes in GRanges object will be used if no sizes there the end
#'        of each chr will be the end last feature on each chr
#' @param stranded if FALSE, all of the returned features will be strandless
#'        (will have "*" in the strand slot)
#' @param keep.strand.prop If TRUE strands will have the same proportion as the
#'        features
#' @param keep.chrom If TRUE, number of features and randomized features for a
#'        chromosome will match. Currently seeting this to FALSE is not supported.
#' @param exclude A GRanges object where no randomized feature should overlap,
#'        can be gaps or unmappable regions in the genome as an example.
#' @param include A GRanges object which defines the boundaries of randomized
#'        features. If not provided the whole genome is used, as defined using the chrom.sizes parameter.
#' @param seed random number generator seed
#' @param nrand number of randomizations (default:1)
#' @usage randomizeFeature(feature,chrom.sizes=NULL,stranded=TRUE,
#'                          keep.strand.prop=TRUE,keep.chrom=TRUE,
#'                          exclude=NULL,include=NULL,seed=NULL)
#' @return returns a GRanges object which is randomized version of the feature, along with a "set" column in the metadata which designates to which iteration of the randomization the range belong.

#' @export
#' @docType methods
#' @rdname randomizeFeature-methods
setGeneric("randomizeFeature", function(feature,chrom.sizes=NULL,stranded=TRUE,
keep.strand.prop=TRUE,keep.chrom=TRUE,
exclude=NULL,include=NULL,seed=NULL,nrand=1)
standardGeneric("randomizeFeature") )

#' @aliases randomizeFeature,GRanges-method
#' @rdname randomizeFeature-methods
setMethod("randomizeFeature", signature(feature = "GRanges"),
			function( feature,
					  chrom.sizes ,
					  stranded ,
					  keep.strand.prop, 
					  keep.chrom,
					  exclude,
					  include,
					  seed,
					  nrand ){

			
			
			# sets the random seed generator if it isn't set
			if(!is.null(seed) & is.numeric(seed)){set.seed(seed)}
			
			#checks whether the chromosome sizes are defined, if not uses the end coordinates of features
			if(is.null(chrom.sizes)){
				cat('Using feature chromosome sizes ...')
				chrom.sizes=list()
				chrs=as.character(unique(seqnames(feature)))
				my.starts=c()
				my.ends=c()
				for( i in 1:length(chrs) ){
					my.set= feature[seqnames(feature)==chrs[i], ]
					chrom.sizes[chrs[i]]  =max( end(my.set) )
				}
				chrom.sizes=unlist(chrom.sizes)
			}

			# if the include ranges are not defined, defines the ranges across the whole chromosome
			if(is.null(include))
				include = GRanges(names(chrom.sizes), IRanges(1, chrom.sizes))
			
			# if the exclude ranges are not defined it defines them
			if(is.null(exclude))
				exclude = GRanges()
			
			# gives a warning if the restricted ranges are too small
			if(sum(width(feature))/(sum(as.numeric(width(include))) - sum(as.numeric(width(exclude)))) < 1e-3)
				warning('The width of the feature sizes is too big compared to the allowed genome')
			
			# checks whether the features are reduced
			if(sum(width(feature)) != sum(width(reduce(feature))))
				stop('Feature needs to be a disjoined set of ranges')
			
			# selects the chromosome set to be used
			chrs = intersect(as.character(seqnames(include)), names(chrom.sizes))
			chrom.sizes = chrom.sizes[chrs]
			exclude = exclude[seqnames(exclude) %in% chrs]
			
			# defines the regions that can be sampled from
			cat('Creating allowed regions...\n')
			r = sapply(chrom.sizes, function(x)Rle(0, x))
			re = sapply(names(r), function(x){a=r[[x]];a[ranges(include[seqnames(include) == x])]= 1;a})
			if(length(exclude) > 0)
				re = sapply(names(r), function(x){a=re[[x]];a[ranges(exclude[seqnames(exclude) == x])]= 0;a})
			names(re) = names(chrom.sizes)
			
			# loops over the chromosome and constructs the defined set of ranges
			cat('Looping over the chromosomes...\n')
			glist = list()
			for(chr in chrs){
			
				cat(chr,'\r')
				# gets the indices of the features on the chromosome
				ind = seqnames(feature) == chr
				
				# gets the permitted ranges
				a = ranges(re[[chr]])[runValue(re[[chr]]) == 1]
				nl = sum(ind)
				total = nrand*nl
				
				# gets the weight for the number of regions per allowed region based on the width
				w = round(total*(width(a)/sum(width(a))))
				i = sample(length(w),1)
				w[i] = w[i] + total - sum(w)
				
				# generates the random regions
				m = mapply(function(n,x,y)runif(n,x,y), w, start(a), end(a), SIMPLIFY=F)
				g = GRanges(chr, IRanges(start = unlist(m), width=rep(width(feature[ind]), times=nrand)))
				values(g)$set = rep(1:nrand, each=nl)
				
				# sets the strand of the features
				if(stranded == TRUE)
					strand(g) = sample(c('+','-'), length(g), replace=T)
				if(keep.strand.prop == TRUE)
					strand(g) = sample(as.character(strand(feature[ind])))
				glist[[chr]] = g
			}
			cat('Returning the final GRanges object...\n')	
			return(unlist(GRangesList(glist)))
})


# ------------------------------------------------------------------------------------------------ #
# calculateOverlapSignificance
#' function that calculates the significance of overlaps of two sets of features using randomization
#'
#' This function calculates the significance of overlaps of two sets of features using randomization. #' It returns a distributon of overlaps of a target set with a given randomized feature set. The 
#' randomization can be constrained by supplied arguments.
#' The function is still in Beta mode - the regions can overlap excluded regions, and the randomized regions are not disjoint. Please take care that the excluded and included regions are not too strict when compared to the total width of the ranges.
#'
#' @param target a GRanges object for which the overlap needs to be calculates
#' @param feature a GRanges object to be randomized
#' @param chrom.sizes sizes of chromosomes as a named vector (names are
#'        chromsomes names and elements of the vectors are lengths). , if not
#'        given sizes in GRanges object will be used if no sizes there the end
#'        of each chr will be the end last feature on each chr
#' @param stranded if FALSE, all of the returned features will be strandless
#'        (will have "*" in the strand slot)
#' @param keep.strand.prop If TRUE strands will have the same proportion as the
#'        features
#' @param keep.chrom If TRUE, number of features and randomized features for a
#'        chromosome will match. Currently seeting this to FALSE is not supported.
#' @param exclude A GRanges object where no randomized feature should overlap,
#'        can be gaps or unmappable regions in the genome as an example.
#' @param include A GRanges object which defines the boundaries of randomized
#'        features
#' @param seed random number generator seed
#' @param nrand number of randomizations (default:1)
#' @usage randomizeFeature(feature,chrom.sizes=NULL,stranded=TRUE,
#'                          keep.strand.prop=TRUE,keep.chrom=TRUE,
#'                          exclude=NULL,include=NULL,seed=NULL)
#' @return returns a GRanges object which is randomized version of the feature
#' @export
#' @docType methods
#' @rdname randomizeFeature-methods
setGeneric("calculateOverlapSignificance", 
				function(target, feature, chrom.sizes=NULL, stranded=TRUE,
						 keep.strand.prop=TRUE, keep.chrom=TRUE,
						 exclude=NULL, include=NULL, seed=NULL, nrand=1)
				standardGeneric("calculateOverlapSignificance") )

setMethod("calculateOverlapSignificance", signature(target="GRanges", feature="GRanges"),
			function( target,
					  feature,
					  chrom.sizes ,
					  stranded ,
					  keep.strand.prop, 
					  keep.chrom,
					  exclude,
					  include,
					  seed,
					  nrand ){
		
		# randomizes the ranges nrand times
		rand.ranges = randomizeFeature(feature=feature, 
										chrom.sizes=chrom.sizes,
										stranded=stranded, 
										keep.strand.prop=keep.strand.prop, 
										keep.chrom=keep.chrom, 
										exclude=exclude, 
										include=include, 
										seed=seed,nrand = nrand)
		
		# converts the overlaps into a data.table object and summarizes the overlaps
		cat('Summarizing the overlaps...\n')
		co = data.table(as.matrix(findOverlaps(rand.ranges, target)))
		co$set = values(rand.ranges)$set[co$queryHits]
		co$queryHits = NULL
		co = unique(co)
		co = co[,length(subjectHits), by=set]
		return(co$V1/length(target))
})
