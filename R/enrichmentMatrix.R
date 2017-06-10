#' Compute an enrichment of IP over control stored in ScoreMatrixControl object 
#' 
#' This is a \code{\link{enrichmentMatrix}} function for a ScoreMatrixControl object that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' @param x the \code{\link{ScoreMatrixControl}} object 
#' @return \code{ScoreMatrix} 
#' @aliases enrichmentMatrix,ScoreMatrixControl-method
#' @rdname enrichmentMatrix-ScoreMatrixControl-method
setGeneric("enrichmentMatrix",
           function(x) 
             standardGeneric("enrichmentMatrix") )


setMethod("enrichmentMatrix", signature("ScoreMatrixControl"),
          function(x) {
            enrichment <- matrix()
            enrichment <- log2((x@.Data + 1) / (x@control + 1))
            enrichment <- as(enrichment,"ScoreMatrix")
            validObject(enrichment)
            return(enrichment) 
          }
)
