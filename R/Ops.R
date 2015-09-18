## Ops methods for genomation classes: ScoreMatrix, ScoreMatrixList

# Ops method for a ScoreMatrix object. It enables to use arithmetic, indicator and logic operations on ScoreMatrix objects.
setMethod("Ops", signature(e1="ScoreMatrix", e2="ScoreMatrix"),
          function(e1, e2) {
            e1@.Data=callGeneric(e1@.Data, e2@.Data)
            validObject(e1)
            return(e1)
          }
)

# Subseting method for a ScoreMatrix object.
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

# Ops method for a ScoreMatrixList object. It enables to use arithmetic, indicator and logic operations on ScoreMatrixList objects.
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
              smllist[[i]] <- callGeneric(e1d[[i]], e2[[i]])
            }
            sml = as(smllist, "ScoreMatrixList")
            validObject(sml)
            return(sml)
          }
)

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

# Subseting method for a ScoreMatrixList object.
setMethod("[",signature(x="ScoreMatrixList", i = "ANY"), 
          function(x,i){
            ScoreMatrixList(x@.Data[i])
          }
)

