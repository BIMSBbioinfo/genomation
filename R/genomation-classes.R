# S4 classes for genomation show and accessor functions

# ---------------------------------------------------------------------------- #
# --------------------- #
# Randomization Classes #
# --------------------- #
#' An S4 class for storing \code{getRandomEnrichment} function results
#'
#' The resulting object stores the results of \code{getRandomEnrichment} function
#'
#' @section Slots:\describe{
#'                  \item{\code{orig.cnt}:}{number of features overlapping with query at \code{getRandomEnrichment} }
#'                  \item{\code{rand.olap.dist}:}{set of number of features overlapping with randomized queries at \code{getRandomEnrichment}}
#'                  \item{\code{log2fc}:}{ log2 fold change calculated by dividing \code{orig.cnt} by mean(\code{rand.olap.dist}) and taking log2 of that result}
#'                  \item{\code{p.value}:}{P-value assuming \code{rand.olap.dist} has a normal distribution and comparing \code{orig.cnt} with that distribution }
#'                  \item{\code{rand.p.value}:}{ p-value from randomization by calculation the proportion of how many times a random number of overlap exceeds the original number of overlap}
#'                 }
#'
#' @name RandomEnrichment-class
#' @rdname RandomEnrichment-class
#' @seealso \code{\link{getRandomEnrichment}}
#' @export
setClass("RandomEnrichment", 
		representation(
			orig.cnt = "numeric", 
			rand.olap.dist = "numeric",
			log2fc="numeric",
			p.value="numeric",
			rand.p.value="numeric" 
))


#' show method for some of the genomation classes
#' @param object object of class RandomEnrichment
#' @rdname show-methods
#' @aliases show,RandomEnrichment-method
setMethod("show", "RandomEnrichment", function(object) {
  message("orig.cnt:",object@orig.cnt,"\n")
  message("log2fc:",object@log2fc,"\n")
  message("p.value:",object@p.value,"\n")
  message("rand.p.value:", object@rand.p.value,"\n")
})





# ---------------------------------------------------------------------------- #
# ------------------- #
# ScoreMatrix Classes #
# ------------------- #

#' An S4 class for storing \code{ScoreMatrix} function results
#'
#' The resulting object is an extension of a \code{matrix} object, and stores values (typically genome-wide scores) for a predefined set of regions
#' Each row on the ScoreMatrix is a predefined region (Ex: CpG islands, promoters) and columns are values across those regions.
#'
#' @name ScoreMatrix-class
#' @rdname ScoreMatrix-class
#' @seealso \code{\link{ScoreMatrix}}
#' @export
setClass("ScoreMatrix",contains = "matrix")


# ---------------------------------------------------------------------------- #
#' An S4 class for storing a set of \code{ScoreMatrixList} 
#'
#' The resulting object is an extension of a \code{list} object, where each element corresponds to a score matrix object
#'
#' @name ScoreMatrixList-class
#' @rdname ScoreMatrixList-class
#' @seealso \code{\link{ScoreMatrixList}}
#' @export
setClass("ScoreMatrixList", 
			contains = "list"
#			validity=.valid.ScoreMatrixList
			)


# ---------------------------------------------------------------------------- #
# ------------------ #
# Annotation Classes #
# ------------------ #

# A set of objects that will hold statistics about feature and annotation overlap
#' An S4 class that information on overlap of target features with annotation features  
#'
#' This object is desgined to hold statistics and information about genomic feature overlaps
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target features with annotation genomic features}
#'
#'                  \item{\code{annotation}}{a named vector of percentages}
#'
#'                  \item{\code{precedence}}{a named vector of percentages}
#'
#'                  \item{\code{num.hierarchica}}{vector}
#'
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#' }
#' @name AnnotationByFeature-class
#' @rdname AnnotationByFeature-class
#' @export
setClass("AnnotationByFeature", 
         representation(members  = "matrix",
                        annotation = "numeric",
                        precedence = "numeric",
                        num.annotation = "numeric",
                        num.precedence = "numeric",
                        no.of.OlapFeat = "numeric",
                        perc.of.OlapFeat = "numeric"))



#' An S4 class that information on overlap of target features with annotation features  
#'
#' This object is desgined to hold statistics and information about genomic feature overlaps
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target features with annotation genomic features}
#'
#'                  \item{\code{annotation}}{a named vector of percentages}
#'
#'                  \item{\code{precedence}}{a named vector of percentages}
#'
#'                  \item{\code{num.hierarchica}}{vector}
#'
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#'
#'                  \item{dist.to.TSS}{a data frame showing distances to TSS and gene/TSS names and strand}
#' }
#' @name AnnotationByGeneParts-class
#' @rdname AnnotationByGeneParts-class
#' @export
setClass("AnnotationByGeneParts", 
         representation(dist.to.TSS = "data.frame"), contains = "AnnotationByFeature")



