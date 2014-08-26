# ---------------------------------------------------------------------------- #
#' get enrichment based on randomized feature overlap
#'
#' This function measures the association between two genomic features by randomizing one feature and counting the overlaps in randomized sets.
#' That is to say, \code{query} feature will be randomly distributed over the genome (constrained by provided options), and the overlap of \code{target} with these randomized features will be measured.
#'
#' @param target a \code{GRanges} object to be overlapped with \code{query}
#' @param query a \code{GRanges} object that will be randomly placed across the genome and overlap of these random regions with \code{target} will be the background distribution of association between \code{target} and \code{query}.
#' @param rand.set instead of randomly placing features in \code{query} one can supply an already shuffled set of \code{query} genomic features.
#' @param randomizations  number of times the features to be shuffled
#' @param ... other parameters to be passed to \code{randomizeFeature} function. These parameters ccontrol how randomization is done.
#' @return returns a \code{RandomEnrichment} object
#' @seealso \code{\link{randomizeFeature}}
#' @examples
#' # data(cage)
#' # data(cpgi)
#'  
#' # enr = getRandomEnrichment(cage, cpgi, randomizations=50)
#' 
#' 
#' @export
#' @docType methods
#' @rdname getRandomEnrichment-methods
setGeneric("getRandomEnrichment", 
           function(target, query, randomizations=1000, rand.set=NULL,...) 
             standardGeneric("getRandomEnrichment") )

#' @aliases getRandomEnrichment,GRanges,GRanges-method
#' @rdname getRandomEnrichment-methods
setMethod("getRandomEnrichment",
			signature(target = "GRanges",query="GRanges"),
			function( target, query, randomizations, rand.set,... ){

				orig.cnt = sum(countOverlaps(target,query) > 0)

				if( is.null(rand.set) ){
				
					rand.olap.dist=numeric(randomizations)
					for(i in 1:randomizations){
# 						message("Iteration number:",i,"\r")
						my.rand = randomizeFeature(query,...)
						my.cnts = sum(countOverlaps(target,my.rand)>0)
						rand.olap.dist[i] = my.cnts
					}
# 					message("\n")
					
				}else if(is(rand.set,"GRangesList")){
					
					randomizations = length(rand.set)
					rand.olap.dist = numeric(randomizations)
					for(i in 1:randomizations){
						message("Iteration number:",i,"\r")
						my.cnts = sum(countOverlaps(target,rand.set[[i]])>0)
						rand.olap.dist[i] = my.cnts
					}
					message("\n")
				}else{
					stop("Wrong 'rand.set' argument supplied, it should be a 'GRangesList'")
				}

				# calculate p-value assuming normal distribution
				p.value=pnorm(orig.cnt, 
                      mean=mean(rand.olap.dist), 
                      sd=sd(rand.olap.dist),
                      lower.tail=FALSE )
				#if(p1>0.5){p.value=1-p1}else{p.value=p1}

				# calculate randomized raw p-value
				rand.p.value = sum(rand.olap.dist>orig.cnt)/length(rand.olap.dist)

				n = new('RandomEnrichment', 
                orig.cnt = orig.cnt,
					      rand.olap.dist = rand.olap.dist,
				        log2fc = log2(orig.cnt/mean(rand.olap.dist)),
					      p.value = p.value,
					      rand.p.value = rand.p.value)
        return(n)
})
