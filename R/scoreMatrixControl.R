# ---------------------------------------------------------------------------- #
#' Make ScoreMatrixControl from ScoreMatrix objects of IP and control samples
#' 
#' The function creates a \code{ScoreMatrixControl} object that stores a \code{ScoreMatrix} containing a signal values of IP sample as well as a \code{ScoreMatrix} object with control sample. 
#' This object is an input in \code{\link{enrichmentMatrix}} function that compute enrichment over IgG or input DNA control.
#'
#' @param IP the \code{\link{ScoreMatrix}} object storing an IP sample
#' @param control the \code{\link{ScoreMatrix}} object storing a control sample
#' @return returns a \code{ScoreMatrixControl} object
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
            m <- as(IP, "ScoreMatrixControl")
            m@control <- control
            validObject(m)
            return(m)
          }
)  