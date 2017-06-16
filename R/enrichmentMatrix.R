#' Compute an enrichment of IP over control both stored in ScoreMatrixControl object
#' 
#' This is a \code{enrichmentMatrix} function for \code{ScoreMatrixControl} object, 
#' that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' 
#' @param x \code{\link{ScoreMatrixControl}}
#' @return \code{\link{ScoreMatrix}}
#' @aliases enrichmentMatrix,ScoreMatrixControl-method
#' @rdname enrichmentMatrix-ScoreMatrixControl-method
#' @seealso \code{\link{ScoreMatrixListControl}}, \code{\link{ScoreMatrixControl}}
#' @examples  
#' #load IP and control BAM files and create ScoreMatrix objects
#' source("http://bioconductor.org/biocLite.R")
#' biocLite('genomationData')
#' bam.file_IP <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescSuz12051317AlnRep1.chr21.bam", package = "genomationData")
#' bam.file_c <- system.file("extdata", 
#' "wgEncodeBroadHistoneH1hescControlStdAlnRep1.chr21.bam", package = "genomation")
#' data(promoters)
#' IP <- ScoreMatrix(target = bam.file_IP, windows = promoters, type = 'bam')
#' control <- ScoreMatrix(target = bam.file_c, windows = promoters, type = 'bam')
#' 
#' # create a ScoreMatrixControl object
#' smControl <- ScoreMatrixControl(IP, control)
#' 
#' # compute an enrichment of IP over control
#' enrichmentMatrix(smControl)
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
#' This is a \code{enrichmentMatrix} function for \code{ScoreMatrixListControl} object, 
#' that enables to normalize ChIP-seq signals with respect to IgG or input DNA control.
#' 
#' @param x \code{\link{ScoreMatrixListControl}} 
#' @return \code{\link{ScoreMatrixList}}
#' @aliases enrichmentMatrix,ScoreMatrixListControl-method
#' @rdname enrichmentMatrix-ScoreMatrixListControl-method
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
#' "wgEncodeBroadHistoneH1hescControlStdAlnRep1.chr21.bam", package = "genomation")
#' control <- ScoreMatrix(target = bam.file_c, windows = promoters, type = 'bam')
#' 
#' # create a ScoreMatrixList object of IP ScoreMatrix object
#' IPs <- list(IP_1,IP_2)
#' sml_IP <- ScoreMatrixList(IPs)
#' 
#' #create a ScoreMatrixListControl object
#' smlControl <- ScoreMatrixListControl(sml_IP, control)
#' 
#' # compute an enrichment of IP over control
#' enrichmentMatrix(smlControl)
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