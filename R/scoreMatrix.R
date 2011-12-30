

#' get the scores of
#'
#' This function measures the association between two genomic features by randomizing one feature and counting the overlaps in randomized sets. 
#' That is to say, \code{query} feature will be randomly distributed over the genome (constrained by provided options), and the overlap of \code{target} with these randomized features will be measured.
#'
#' @param target a \code{RleList} or a \code{modRleList} object to be overlapped with ranges in \code{windows}
#' @param windows a \code{GRanges} object that will be randomly placed across the genome and overlap of these random regions with \code{target} will be the background distribution of association between \code{target} and \code{query}.
#'
#' @usage scoreMatrix(target,windows)
#' @return returns a \code{scoreMatrix} object
#' @export
#' @docType methods
#' @rdname scoreMatrix-methods           
setGeneric("scoreMatrix",function(target,windows,...) standardGeneric("scoreMatrix") )

#' @alias scoreMatrix,GRanges,RleList-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("GRanges","RleList"),
          function(target,windows,...){
            
            #check if all windows are equal length
            
            #get views 
            
            # make a matrix from those views
            
})

#' @alias scoreMatrix,GRanges,RleList-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("GRanges","GRanges"),
          function(target,windows){
            
            #check if all windows are equal length
            
            #get views 
            
            # make a matrix from those views
            
})