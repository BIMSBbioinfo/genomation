# ---------------------------------------------------------------------------- #
#' @title \code{Ops} method for a \code{ScoreMatrix} object. It enables to use arithmetic, indicator and logic operations on \code{ScoreMatrix} objects.
#' @param e1 the \code{\link{ScoreMatrix}} object or numeric value
#' @param e2 the \code{\link{ScoreMatrix}} object or numeric value
#' @return \code{ScoreMatrix} 
#' @aliases Ops,ScoreMatrix-method
#' @rdname ScoreMatrix-methods
setMethod("Ops", signature(e1="ScoreMatrix", e2="ScoreMatrix"),
          function(e1, e2) {
            e1@.Data=callGeneric(e1@.Data, e2@.Data)
            validObject(e1)
            return(e1)
          }
)