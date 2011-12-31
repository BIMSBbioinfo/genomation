

#' get scores overlapping with windows in a scoreMatrix object
#'
#' A scoreMatrix object can be used to draw average profiles or heatmap of read coverage or wig track-like data.
#' \code{windows} can be a predefined region around transcription start sites or other regions of interest that have equal lengths
#'
#' @param target a \code{RleList} or a \code{modRleList} object to be overlapped with ranges in \code{windows}
#' @param windows a \code{GRanges} object that will be randomly placed across the genome and overlap of these random regions with \code{target} will be the background distribution of association between \code{target} and \code{query}.
#' @param ... parameters to be passed to \code{modCoverage} function
#'
#' @usage scoreMatrix(target,windows)
#' @return returns a \code{scoreMatrix} object
#' @export
#' @docType methods
#' @rdname scoreMatrix-methods           
setGeneric("scoreMatrix",function(target,windows,...) standardGeneric("scoreMatrix") )

#' @aliases scoreMatrix,GRanges,RleList-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("GRanges","RleList"),
          function(target,windows){
            
            #check if all windows are equal length
            
            #get views 
            
            # make a matrix from those views
            
            # if the target is modRleList do appropriate calculations to get the score and put NAs in cells that have no value
            
})

#' @aliases scoreMatrix,GRanges,GRanges,ANY-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("GRanges","GRanges"),
          function(target,windows,...){
            
            
            #check if all windows are equal length
            
            #make coverage vector (modRleList) from target
            
            #get views 
            
            # make a matrix from those views
            
            # if the target is modRleList do appropriate calculations to get the score and put NAs in cells that have no value

            
})



