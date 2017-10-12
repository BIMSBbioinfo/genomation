#' Compute an enrichment of IP over control both stored in ScoreMatrix objects
#' 
#' This is an \code{enrichmentMatrix} function for \code{ScoreMatrix} objects, 
#' that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' @param IP \code{\link{ScoreMatrix}} object storing an IP sample
#' @param control \code{\link{ScoreMatrix}} object storing a control sample
#' @note 
#' The function computes an enrichment of IP over control as follow:
#' Suppose both IP and control are ScoreMatrix objects that have same dimensions. 
#' Then, the enrichment is calculated usign a formula: log2((IP + 1) / (control + 1)).
#' @return \code{ScoreMatrix} object
#' @usage  \\S4method{enrichmentMatrix}{ScoreMatrix,ScoreMatrix}(IP, control)
#' @aliases enrichmentMatrix,ScoreMatrix,ScoreMatrix-method
#' @rdname enrichmentMatrix-ScoreMatrix-method
#' @seealso \code{\link{ScoreMatrix}}
#' @docType methods
#' @examples  
#' #load IP and control BAM files and create ScoreMatrix objects
#' source("http://bioconductor.org/biocLite.R")
#' biocLite('genomationData')
#' bam.file_IP <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescSuz12051317AlnRep1.chr21.bam", package = "genomationData")
#' bam.file_c <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescCtcfStdAlnRep1.chr21.bam", package = "genomationData")
#' data(promoters)
#' IP <- ScoreMatrix(target = bam.file_IP, windows = promoters, type = 'bam')
#' control <- ScoreMatrix(target = bam.file_c, windows = promoters, type = 'bam')
#' 
#' # compute an enrichment of IP over control
#' enrichmentMatrix(IP, control)
#' 
setGeneric("enrichmentMatrix",
           function(IP, control)
             standardGeneric("enrichmentMatrix") )

setMethod("enrichmentMatrix", signature("ScoreMatrix", "ScoreMatrix"),
          function(IP, control) {
            enrichment <- matrix()
            enrichment <- log2((IP + 1) / (control + 1))
            enrichment <- as(enrichment, "ScoreMatrix")
            validObject(enrichment)
            return(enrichment) 
          }
)

# ---------------------------------------------------------------------------- #
#' Compute an enrichment of IP over control both stored in ScoreMatrixList objects  
#' 
#' This is an \code{enrichmentMatrix} function for \code{ScoreMatrixList} objects, 
#' that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' 
#' @param IP \code{\link{ScoreMatrixList}} object storing IP samples
#' @param control \code{\link{ScoreMatrixList}} storing control samples
#' @note 
#' The function computes an enrichment of IP over control as follow:
#' Suppose both IP and control are ScoreMatrix objects that have same dimensions. 
#' Then, the enrichment is calculated usign a formula: log2((IP + 1) / (control + 1)).
#' @return \code{ScoreMatrixList} object
#' @usage  \\S4method{enrichmentMatrix}{ScoreMatrixList,ScoreMatrixList}(IP, control)
#' @aliases enrichmentMatrix,ScoreMatrixList,ScoreMatrixList-method
#' @rdname enrichmentMatrix-ScoreMatrixList-method
#' @seealso \code{\link{ScoreMatrixList}}
#' @docType methods
#' @examples
#' #load IP and control BAM files and create ScoreMatrix objects
#' source("http://bioconductor.org/biocLite.R")
#' biocLite('genomationData')
#' data(promoters)
#' bam.file_IP_1 <- system.file("extdata", 
#' "wgEncodeSydhTfbsH1hescZnf143IggrabAlnRep1.chr21.bam", package = "genomationData")
#' IP_1 <- ScoreMatrix(target = bam.file_IP_1, windows = promoters, type = 'bam')
#'
#' bam.file_IP_2 <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescSuz12051317AlnRep1.chr21.bam", package = "genomationData")
#' IP_2 <- ScoreMatrix(target=bam.file_IP_2, windows = promoters, type = 'bam')
#' 
#'bam.file_c <- system.file("extdata", 
#'"wgEncodeBroadHistoneH1hescCtcfStdAlnRep1.chr21.bam", package = "genomationData")
#' control <- ScoreMatrix(target = bam.file_c, windows = promoters, type = 'bam')
#' 
#' # create a ScoreMatrixList object storing IP ScoreMatrix objects
#' sml_IP <- ScoreMatrixList(list(IP1 = IP_1, IP2 = IP_2))
#' 
#' # create a ScoreMatrixList object storing control ScoreMatrix objects
#' sml_control <- ScoreMatrixList(list(c1 = control, c2 = control))
#' 
#' # compute an enrichment of IP over control
#' enrichmentMatrix(sml_IP, sml_control)
#' 
setMethod("enrichmentMatrix", signature("ScoreMatrixList","ScoreMatrixList"),
          function(IP, control) {
            
            # checks whether the arguments consist of the same number of elements 
            if(length(IP) != length(control))
              stop('IP argument and control argument should have the same length')
            
            # checks whether all corresponding elements of ScoreMatrixList objects have the same size 
            check.size <- function (i, IP, control){
              if(dim(IP[[i]])[1] != dim(control[[i]])[1] || dim(IP[[i]])[2] != dim(control[[i]])[2])
                stop( "element:", i, "\t IP argument and control argument do not have the same size")
            }
            
            n <- mclapply(1:length(IP), check.size, IP = IP, control = control)
            
            calc.enrichment<- function(i, IP, control){
              enrichment <- matrix()
              enrichment <- log2((IP[[i]] + 1) / (control[[i]] + 1))
              enrichment <- as(enrichment, "ScoreMatrix")
            }
            
            smlc <- list()
            smlc <- mclapply(1:length(IP), calc.enrichment, IP=IP, control=control)
            smlc <- as(smlc, "ScoreMatrixList")
            names(smlc) <- names(IP)
            
            validObject(smlc)
            return(smlc) 
          }
)

# -------------------------------------------------------------------------------------------------------- #
#' Compute an enrichment of IP (stored in ScoreMatrixList object) over control (stored in ScoreMatrix object)   
#' 
#' This is an \code{enrichmentMatrix} function for IP \code{ScoreMatrixList} object and control \code{ScoreMatrix} object, 
#' that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' 
#' @param IP \code{\link{ScoreMatrixList}} object storing IP samples
#' @param control \code{\link{ScoreMatrix}} storing control sample
#' @note 
#' The function computes an enrichment of IP over control as follow:
#' Suppose both IP and control are ScoreMatrix objects that have same dimensions. 
#' Then, the enrichment is calculated usign a formula: log2((IP + 1) / (control + 1)).
#' @return \code{ScoreMatrixList} object
#' @usage  \\S4method{enrichmentMatrix}{ScoreMatrixList,ScoreMatrix}(IP, control)
#' @aliases enrichmentMatrix,ScoreMatrixList,ScoreMatrix-method
#' @rdname enrichmentMatrix-ScoreMatrixList-ScoreMatrix-method
#' @seealso \code{\link{ScoreMatrixList}}, \code{\link{ScoreMatrix}}
#' @docType methods
#' @examples
#' #load IP and control BAM files and create ScoreMatrix objects
#' source("http://bioconductor.org/biocLite.R")
#' biocLite('genomationData')
#' data(promoters)
#' bam.file_IP_1 <- system.file("extdata", 
#' "wgEncodeSydhTfbsH1hescZnf143IggrabAlnRep1.chr21.bam", package = "genomationData")
#' IP_1 <- ScoreMatrix(target = bam.file_IP_1, windows = promoters, type = 'bam')
#'
#' bam.file_IP_2 <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescSuz12051317AlnRep1.chr21.bam", package = "genomationData")
#' IP_2 <- ScoreMatrix(target=bam.file_IP_2, windows = promoters, type = 'bam')
#' 
#' bam.file_c <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescCtcfStdAlnRep1.chr21.bam", package = "genomationData")
#' control <- ScoreMatrix(target = bam.file_c, windows = promoters, type = 'bam')
#' 
#' # create a ScoreMatrixList object storing IP ScoreMatrix objects
#' sml_IP <- ScoreMatrixList(list(IP1 = IP_1, IP2 = IP_2))
#'  
#' # compute an enrichment of IP over control
#' enrichmentMatrix(sml_IP, control)
#' 
setMethod("enrichmentMatrix", signature("ScoreMatrixList","ScoreMatrix"),
          function(IP, control) {
            
            # checks whether each element of IP ScoreMatrixList has the same size like control ScoreMatrix 
            check.size <- function (i, IP, control){
              if(dim(IP[[i]])[1] != dim(control)[1] || dim(IP[[i]])[2] != dim(control)[2])
                stop( "element:",i,"\t IP argument and control argument do not have the same size")
            }
            
            n <- mclapply(1:length(IP), check.size, IP = IP, control = control)
            
            calc.enrichment<- function(i, IP, control){
              enrichment <- matrix()
              enrichment <- log2((IP[[i]] + 1) / (control + 1))
              enrichment <- as(enrichment, "ScoreMatrix")
            }
            
            smlc <- list()
            smlc <- mclapply(1:length(IP), calc.enrichment, IP = IP, control = control)
            smlc <- as(smlc, "ScoreMatrixList")
            names(smlc) <- names(IP)
            
            validObject(smlc)
            return(smlc) 
          }
)