# S4 classes for genomation show and accessor functions

# ------------------------------------------------------------------------------ #
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
#' @name randomEnrichment-class
#' @rdname randomEnrichment-class
#' @seealso \code{\link{getRandomEnrichment}}
#' @export
setClass("randomEnrichment", 
		representation(
			orig.cnt = "numeric", 
			rand.olap.dist = "numeric",
			log2fc="numeric",
			p.value="numeric",
			rand.p.value="numeric" 
))


#' show method for some of the genomation classes
#' @rdname show-methods
#' @aliases show,randomEnrichment-method
setMethod("show", "randomEnrichment", function(object) {
  cat("orig.cnt:",object@orig.cnt,"\n")
  cat("log2fc:",object@log2fc,"\n")
  cat("p.value:",object@p.value,"\n")
  cat("rand.p.value:", object@rand.p.value,"\n")
})





# ------------------------------------------------------------------------------ #
#' An S4 class for storing \code{scoreMatrix} function results
#'
#' The resulting object is an extension of a \code{matrix} object, and stores values (typically genome-wide scores) for a predefined set of regions
#' Each row on the scoreMatrix is a predefined region (Ex: CpG islands, promoters) and columns are values across those regions.
#'
#' @name scoreMatrix-class
#' @rdname scoreMatrix-class
#' @seealso \code{\link{scoreMatrix-methods}}
#' @export
setClass("scoreMatrix",contains = "matrix")


# ------------------------------------------------------------------------------ #
#' An S4 class for storing a set of \code{ScoreMatrixList} 
#'
#' The resulting object is an extension of a \code{list} object, where each element corresponds to a score matrix object
#'
#' @name ScoreMatrixList-class
#' @rdname ScoreMatrixList-class
#' @seealso \code{\link{ScoreMatrixList-methods}}
#' @export
setClass("ScoreMatrixList", 
			contains = "list"
#			validity=.valid.ScoreMatrixList
			)