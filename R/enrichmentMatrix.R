#' Compute an enrichment of IP over control both stored in ScoreMatrixControl object
#' 
#' This is a \code{\link{enrichmentMatrix}} function for ScoreMatrixControl object, that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' 
#' @param x \code{\link{ScoreMatrixControl}}
#' @return \code{\link{ScoreMatrix}}
#' @aliases enrichmentMatrix,ScoreMatrixControl-method
#' @rdname enrichmentMatrix-ScoreMatrixControl-method
#' 
setGeneric("enrichmentMatrix",
           function(x)
             standardGeneric("enrichmentMatrix") )

setMethod("enrichmentMatrix", signature("ScoreMatrixControl"),
          function(x) {
            enrichment <- matrix()
            enrichment <- log2((x@.Data + 1) / (x@control + 1))
            enrichment <- as(enrichment, "ScoreMatrix")
            validObject(enrichment)
            return(enrichment) 
          }
)

# ---------------------------------------------------------------------------- #
#' Compute an enrichment of IP over control both stored in ScoreMatrixListControl object  
#' 
#' This is a \code{\link{enrichmentMatrix}} function for ScoreMatrixListControl object, that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' 
#' @param x \code{\link{ScoreMatrixListControl}} 
#' @return \code{\link{ScoreMatrixList}}
#' @aliases enrichmentMatrix,ScoreMatrixListControl-method
#' @rdname enrichmentMatrix-ScoreMatrixListControl-method
#' 
setMethod("enrichmentMatrix", signature("ScoreMatrixListControl"),
          function(x) {
            calc.enrichment<- function(i, x){
              enrichment <- matrix()
              enrichment <- log2((x[[i]]@.Data + 1) / (x[[i]]@control + 1))
              enrichment <- as(enrichment, "ScoreMatrix")
            }
            
            smlc <- list()
            smlc <- mclapply(1:length(x), calc.enrichment, x)
            smlc <- as(smlc, "ScoreMatrixList")
            names(smlc) <- names(x)
            
            validObject(smlc)
            return(smlc) 
          }
)