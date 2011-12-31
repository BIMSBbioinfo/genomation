#' get scores overlapping with windows in a scoreMatrix object
#'
#' A scoreMatrix object can be used to draw average profiles or heatmap of read coverage or wig track-like data.
#' \code{windows} can be a predefined region such as CpG islands or gene bodies that are not necessarily equi-width.
#' Each window will be chopped to equal number of bins based on \code{bin.num} option.
#'
#' @param target a \code{RleList} or a \code{modRleList} object to be overlapped with ranges in \code{windows}
#' @param windows a \code{GRanges} object that will be randomly placed across the genome and overlap of these random regions with \code{target} will be the background distribution of association between \code{target} and \code{query}.
#' @param bin.num A single \code{integer} value denoting how many bins there should be for each window
#' @param ... parameters to be passed to \code{modCoverage} function
#'
#' @usage scoreMatrix(target,windows)
#' @return returns a \code{scoreMatrix} object
#' @export
#' @docType methods
#' @rdname scoreMatrixNorm-methods           
setGeneric("scoreMatrixNorm",function(target,windows,percentile,...) standardGeneric("scoreMatrixNorm") )

#' @alias scoreMatrixNorm,GRanges,RleList-method
#' @rdname scoreMatrixNorm-methods
setMethod("scoreMatrixNorm",signature("GRanges","RleList"),
          function(target,windows,percentile){
            
            #get bins from the windows
            
            #get views for bins and merge adjacent bins
            
            # make a matrix from those views
            
            # if the target is modRleList do appropriate calculations to get the score and put NAs in cells that have no value
            
})

#' @alias  scoreMatrixNorm,GRanges,RleList,ANY-method
#' @rdname scoreMatrixNorm-methods
setMethod("scoreMatrixNorm",signature("GRanges","GRanges"),
          function(target,windows,percentile,...){
            
            
            #get bins from the windows
            
            #make coverage vector (modRleList) from target
            
            #get views for bins and merge adjacent bins
            
            # make a matrix from those views
            
            # if the target is modRleList do appropriate calculations to get the score and put NAs in cells that have no value

            
})
