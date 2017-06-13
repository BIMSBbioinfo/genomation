#' Make ScoreMatrixListControl from ScoreMatrixList objects of IP and control samples
#' 
#' The function creates a \code{ScoreMatrixListControl} object that stores \code{ScoreMatrixControl} objects.
#' This object is an input in \code{\link{enrichmentMatrix}} function that compute enrichment over IgG or input DNA control.
#'
#' @param IP \code{\link{ScoreMatrixList}} object storing IP samples
#' @param control \code{\link{ScoreMatrixList}} storing control samples
#' @return \code{ScoreMatrixListControl} object
#' @aliases ScoreMatrixListControl,ScoreMatrixList,ScoreMatrixList-method
#' @rdname ScoreMatrixListControl-ScoreMatrixList-ScoreMatrixList-method
#' @export
#' @docType methods

setGeneric("ScoreMatrixListControl", 
           function(IP, control){
             standardGeneric("ScoreMatrixListControl" ) 
           })

setMethod("ScoreMatrixListControl", signature("ScoreMatrixList","ScoreMatrixList"),
          function(IP, control){
            
            # checks whether the arguments consist of the same number of elements 
            if(length(IP) != length(control))
              stop('IP argument and control argument should have the same length')
            
            # checks whether all corresponding elements of ScoreMatrixList objects have the same size 
            check.size <- function (i, IP, control){
              if(dim(IP[[i]])[1] != dim(control[[i]])[1] || dim(IP[[i]])[2] != dim(control[[i]])[2])
                stop( "element:",i,"\t IP argument and control argument do not have the same size")
            }
            
            n <- mclapply(1:length(IP), check.size, IP=IP, control=control)
            
            
            #create a ScoreMatrixListControl object
            
            create.ScoreMatrixControl<- function(i, IP, control){
              ScoreMatrixControl(IP[[i]], control[[i]])
            }
            
            smlc <- list()
            smlc <- mclapply(1:length(IP), create.ScoreMatrixControl, IP=IP, control=control)
            names(smlc) <- names(IP)
            smlc <- as(smlc, "ScoreMatrixListControl")
            
            
            validObject(smlc)
            return(smlc)
          }
)  

# ---------------------------------------------------------------------------- #
#' Make ScoreMatrixListControl from ScoreMatrixList object of IP samples and ScoreMatrix object of control sample
#' 
#' The function creates a \code{ScoreMatrixListControl} object that stores \code{ScoreMatrixControl} objects.
#' This object is an input in \code{\link{enrichmentMatrix}} function that compute enrichment over IgG or input DNA control.
#'
#' @param IP \code{\link{ScoreMatrixList}} object storing an IP samples
#' @param control \code{\link{ScoreMatrix}} storing control sample
#' @return \code{ScoreMatrixListControl} object
#' @aliases ScoreMatrixListControl,ScoreMatrixList,ScoreMatrix-method
#' @rdname ScoreMatrixListControl-ScoreMatrixList-ScoreMatrix-method
#' @export
#' @docType methods
setMethod("ScoreMatrixListControl", signature("ScoreMatrixList","ScoreMatrix"),
          function(IP, control){
            
            # checks whether each element of IP ScoreMatrixList has the same size like control ScoreMatrix 
            check.size <- function (i, IP, control){
              if(dim(IP[[i]])[1] != dim(control)[1] || dim(IP[[i]])[2] != dim(control)[2])
                stop( "element:1",i,"\t IP argument and control argument do not have the same size")
            }
            n <- mclapply(1:length(IP), check.size, IP=IP, control=control)
            
            #create a ScoreMatrixListControl object
            
            create.ScoreMatrixControl<- function(i, IP, control){
              ScoreMatrixControl(IP[[i]],control)
            }
            smlc <- list()
            smlc <- mclapply(1:length(IP), create.ScoreMatrixControl, IP=IP, control=control)
            smlc <- as(smlc, "ScoreMatrixListControl")
            names(smlc) <- names(IP)
            
            validObject(smlc)
            return(smlc)
          }
)  

# ---------------------------------------------------------------------------- #
# show Methods
#' @rdname show-methods
#' @return Shows the number of stored ScoreMatrixControl objects and their sizes
#' 
setMethod("show", "ScoreMatrixListControl",
          function(object){
            dims <- lapply(object, dim)
            n <- length(object)
       
            message('ScoreMatrixListControl of length:', n, '\n')
            
            for(i in 1:n){
              s <- sprintf("%d%s %s",i,".","scoreMatrixControl that consists of:")
              message(s)
              
              dimsD = dim(object[[i]]@.Data)
              dimsC = dim(object[[i]]@control)
              s2=sprintf("      %s %d %d", "IP matrix with dims:", dimsD[1], dimsD[2])
              s3=sprintf("      %s %d %d", "control scoreMatrix with dims:", dimsC[1], dimsC[2])
              message(s2,'\n', s3)
            }
          }
)