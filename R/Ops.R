# ---------------------------------------------------------------------------- #
#' Arithmetic method for ScoreMatrix and ScoreMatrixList classes
#' @title Ops method for a ScoreMatrix object. It enables to use arithmetic, indicator and logic operations on ScoreMatrix objects.
#' @param e1 the \code{\link{ScoreMatrix}} object or numeric value
#' @param e2 the \code{\link{ScoreMatrix}} object or numeric value
#' @return \code{ScoreMatrix} 
#' @aliases Ops,ScoreMatrix,ScoreMatrix-method
setMethod("Ops", signature(e1="ScoreMatrix", e2="ScoreMatrix"),
          function(e1, e2) {
            e1@.Data=callGeneric(e1@.Data, e2@.Data)
            validObject(e1)
            return(e1)
          }
)

#' Extract method for a ScoreMatrix object. 
#' 
# @aliases [,ScoreMatrix-method
#' @param x the \code{\link{ScoreMatrix}} object
#' @param i numeric value
#' @param j numeric value
#' @aliases extract,ScoreMatrix,ANY-method
setMethod("[", signature(x="ScoreMatrix", i = "ANY", j="ANY"),  
          function(x,i,j){
            if(missing(j)){
              res=new("ScoreMatrix",x@.Data[i,])
            }
            else if(missing(i)){
              res=new("ScoreMatrix",x@.Data[,j])  
            }
            else{
              res=new("ScoreMatrix",x@.Data[i,j])  
            }
            
            res
          }
)

# ---------------------------------------------------------------------------- #
#' Arithmetic methods for ScoreMatrixList
#' @title Ops method for a ScoreMatrixList object. It enables to use arithmetic, indicator and logic operations on ScoreMatrixList objects.
#' @param e1 the \code{\link{ScoreMatrixList}} object 
#' @param e2 the \code{\link{ScoreMatrixList}} object
#' @return \code{ScoreMatrixList} 
#' @aliases Ops,ScoreMatrixList,ScoreMatrixList-method
setMethod("Ops", signature(e1="ScoreMatrixList", e2="ScoreMatrixList"),
          function(e1, e2) {
            e1d = e1@.Data
            e2d = e2@.Data
            if( length(e1d) != length(e2d) ){
              stop("ScoreMatrixList objects must have the same length")
            }
            #sml = lapply(1:length(e1d), function(i) callGeneric(e1d[[i]],e2d[[i]]))
            # it doesnt work: Error in get(as.character(call[[1L]]), envir = methodEnv) : 
            # object 'FUN' not found 
            smllist = list()
            for(i in 1:length(e1d)){
              smllist[[i]] <- callGeneric(e1d[[i]], e2d[[i]])
            }
            sml = as(smllist, "ScoreMatrixList")
            validObject(sml)
            return(sml)
          }
)

# ---------------------------------------------------------------------------- #
#' Arithmetic method for ScoreMatrixList
#' @title Ops method for a ScoreMatrixList object. It enables to use arithmetic, indicator and logic operations on ScoreMatrixList objects.
#' @param e1 the \code{\link{ScoreMatrixList}} object 
#' @param e2 the numeric value
#' @return \code{ScoreMatrixList} 
#' @aliases Ops,ScoreMatrixList,numeric-method
setMethod("Ops", signature(e1="ScoreMatrixList", e2="numeric"),
          function(e1, e2) {
            e1d <- e1@.Data #list
            smllist = list()
            for(i in 1:length(e1d)){
              smllist[[i]] <- callGeneric(e1d[[i]], e2)
            }
            sml = as(smllist, "ScoreMatrixList")
            validObject(sml)
            return(sml)
          }
)

# ---------------------------------------------------------------------------- #
#' Arithmetic method for ScoreMatrixList
#' @title Ops method for a ScoreMatrixList object. It enables to use arithmetic, indicator and logic operations on ScoreMatrixList objects.
#' @param e1 the numeric value
#' @param e2 the \code{\link{ScoreMatrixList}} object
#' @return \code{ScoreMatrixList} 
#' @aliases Ops,numeric,ScoreMatrixList-method
setMethod("Ops", signature(e1="numeric", e2="ScoreMatrixList"),
          function(e1, e2) {
            e2d <- e2@.Data #list
            smllist = list()
            for(i in 1:length(e2d)){
              smllist[[i]] <- callGeneric(e1, e2d[[i]])
            }
            sml = as(smllist, "ScoreMatrixList")
            validObject(sml)
            return(sml)
          }
)

#' Extract method for a ScoreMatrixList object. 
#' 
# @aliases [,ScoreMatrixList-method
#' @param x the \code{\link{ScoreMatrixList}} object
#' @param i numeric value
#' @aliases extract,ScoreMatrixList,ANY-method
setMethod("[",signature(x="ScoreMatrixList", i = "ANY"), 
          function(x,i){
            tmp=sml@.Data[i]
            names(tmp)=names(x)[i]
            ScoreMatrixList( tmp)
          }
)

