# ---------------------------------------------------------------------------- #
#' Make ScoreMatrixControl from ScoreMatrix objects of IP and control samples
#' 
#' The function creates a \code{ScoreMatrixControl} object that stores a \code{ScoreMatrix} containing a signal values of IP sample as well as a \code{ScoreMatrix} object with control sample. 
#' This object is an input in \code{\link{enrichmentMatrix}} function that compute enrichment over IgG or input DNA control.
#'
#' @param IP \code{\link{ScoreMatrix}} object storing an IP sample
#' @param control \code{\link{ScoreMatrix}} object storing a control sample
#' @return \code{ScoreMatrixControl} object
#' 
#' @aliases ScoreMatrixControl,ScoreMatrix,ScoreMatrix-method
#' @rdname ScoreMatrixControl-ScoreMatrix-ScoreMatrix-method
#' @export
#' @docType methods
setGeneric("ScoreMatrixControl", 
           function(IP, control){
           standardGeneric("ScoreMatrixControl" ) 
})


setMethod("ScoreMatrixControl", signature("ScoreMatrix", "ScoreMatrix"),
          function(IP, control) {
            
            if(dim(IP)[1] != dim(control)[1] || dim(IP)[2] != dim(control)[2]) 
              stop('IP argument and control argument should have the same size')
            
            m <- as(IP,"ScoreMatrixControl")
            m@control <- control
            validObject(m)
            return(m)
          }
)  

# ---------------------------------------------------------------------------- #
# show Methods
#' @rdname show-methods
#' @return Shows the sizes of both IP and control matrices stored in a ScoreMatrixControl object
setMethod("show", "ScoreMatrixControl",
          function(object){
            message('ScoreMatrixControl that consists of:')
            dims <- dim(object@.Data)
            dimsC <- dim(object@control)
            s <- sprintf("  %s %d %d", "IP matrix with dims:", dims[1], dims[2])
            s2 <- sprintf("  %s %d %d", "control scoreMatrix with dims:", dimsC[1], dimsC[2])
            message(s,'\n', s2)  
          })
