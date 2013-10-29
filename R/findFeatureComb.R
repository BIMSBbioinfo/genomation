# ------------------------------------------------------------------------------------ #
#' Provided a GRangesList, finds the combinations of sets of ranges. It is mostly used to look at the combinatorics of transcription factor binding. The function works by, firstly, constructing a union of all ranges in the list, which are then designated by the combinatorics of overlap with the original sets. A caveat of this approach is that the number of possible combinations increases exponentially, so we would advise you to use it with up to 6 data sets. If you wish to take a look at a greater number of factors, methods like self organizing maps or ChromHMM might be more appropriate.

#' @param gl a \code{GRangesList} object, containing ranges for which represent regions enriched for transcription factor binding
#' @param width \code(integer) is the requested width of each enriched region. If 0 the ranges are not resized, if a positive integer, the width of all ranges is set to that number. Ranges are resized relative to the center of original ranges.
#' @param use.names a boolean which tells the function whether to return the resulting ranges with a numeric vector which designates each class (the default), or to construct the names of each class using the names from the GRangesList
#' @param collapse.char a character which will be used to separate the class names if use.names=TRUE. The default is ':'

#' @usage findTFomb(gl, width, use.names=FALSE, collapse.char=':')
#' @return a GRanges object

#' @docType methods
#' @rdname findFeatureComb-methods
#' @export
setGeneric("findFeatureComb", function(gl, width=0, use.names=FALSE, collapse.char=':') standardGeneric("findFeatureComb") )

#' @aliases findFeatureComb
#' @rdname findFeatureComb-methods
setMethod("findFeatureComb", signature("GRangesList"),
          function(gl, width, use.names, collapse.char){
          
              # returns a vector of ones if there is only one ranges set
              if(length(gl) == 1)
                return(rep(1, length(gl[[1]])))
              
              if(width < 0)
                stop('width needs to be a positive integer')
              
              # resizes the ranges before doing the overlaps
              if(width > 0)
                gl = endoapply(gl, function(x)resize(x, width=width, fix='center'))
                
              # makes a union of ranges and calculates the combinatorics
              r = reduce(unlist(gl))
              do = data.frame(lapply(gl, function(x)as.numeric(countOverlaps(r, x) > 0)))
              
              if(use.names == FALSE){
                values(r)$class = as.numeric(factor(rowSums(t(t(do)*(2^(0:(ncol(do)-1)))))))
                return(r)
              
              }
              if(use.names == TRUE){
                # finds the combinations of ranges for each interval in the reduced ranges
                do.comb = as.numeric(factor(rowSums(t(t(do)*(2^(0:(ncol(do)-1)))))))
                cnames = colnames(do)
                
                # selects only unique representatives of combinations
                do = do[!duplicated(do.comb),]
                # pastes names of individual sets
                cnames = apply(do, 1, function(x)paste(cnames[x == 1]), collapse=collapse.char)
                
                
                r$class = cnames[match(do.comb, unique(do.comb))]
                return(r)
              }
          }
)