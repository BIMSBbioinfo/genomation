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
#' 
#' @section Constructors:
#' see \code{\link{ScoreMatrix}}
#' 
#' @section Coercion:
#' as(from, "matrix"): Creates a matrix from \code{ScoreMatrix} object. You can also use S3Part() function to extract the matrix from \code{ScoreMatrix} object.
#' 
#' @section Subsetting:
#' In the code snippets below, x is a ScoreMatrix object.
#' \code{'x[i,j]'}: Get or set elements from row \code{i} and column \code{j} and return a subset ScoreMatrix object.
#' 
#' @seealso \code{\link{ScoreMatrix}}
#' @export
setClass("ScoreMatrix",contains = "matrix")

# ---------------------------------------------------------------------------- #
#' An S4 class for storing \code{ScoreMatrixControl} function results
#'
#' The resulting object is an extension of a \code{ScoreMatrix} object, that stores ScoreMatrix of IP sample as well as ScoreMatrix containing IgG or input DNA control sample. 
#' 
#' @name ScoreMatrixControl-class
#' @rdname ScoreMatrixControl-class
#'  
#' @section Constructors:
#' see \code{\link{ScoreMatrixControl}}
#' 
#' @section Coercion:
#' as(from, "ScoreMatrixControl"): Creates a \code{ScoreMatrixControl} object from a \code{\link{ScoreMatrix}} or \code{\link{ScoreMatrixBin}} objects. 
#' 
#' @section Subsetting:
#' In the code snippets below, x is a ScoreMatrixControl object.
#'  
#' \code{'x[i,j]'}: Get or set elements from row \code{i} and column \code{j} and return a subset of IP ScoreMatrix object.
#' 
#' \code{'x@control[i,j]'}: Get or set elements from row \code{i} and column \code{j} and return a subset of control ScoreMatrix object.
#' 
#' @slot control A \code{\link{ScoreMatrix}} object storing a control sample.
#' 
#' @seealso \code{\link{ScoreMatrixListControl}}, \code{\link{ScoreMatrix}}
#' @export
setClass("ScoreMatrixControl",
         slots = c(control = "ScoreMatrix"),
         contains = "ScoreMatrix"
)

# ---------------------------------------------------------------------------- #
#' An S4 class for storing a set of \code{ScoreMatrixList} 
#'
#' The resulting object is an extension of a \code{list} object, where each element corresponds to a score matrix object
#' 
#' @name ScoreMatrixList-class
#' @rdname ScoreMatrixList-class
#'  
#' @section Constructors:
#' see \code{\link{ScoreMatrixList}}
#' 
#' @section Coercion:
#' as(from, "ScoreMatrixList"): Creates a \code{ScoreMatrixList} object from a list containing \code{\link{ScoreMatrix}} or \code{\link{ScoreMatrixBin}} objects. 
#' 
#' @section Subsetting:
#' In the code snippets below, x is a ScoreMatrixList object.
#' 
#' \code{x[[i]]},\code{x[[i]]}: Get or set elements \code{i}, where \code{i} is a numeric or character vector of length 1.
#' 
#' \code{x$name}, \code{x$name}: value: Get or set element \code{name}, where \code{name} is a name or character vector of length 1.
#' 
#' @seealso \code{\link{ScoreMatrixList}}
#' @export
setClass("ScoreMatrixList",
            slots = c(names = "character"),
            contains = "list"
#            validity=.valid.ScoreMatrixList
            )

# ---------------------------------------------------------------------------- #
#' An S4 class for storing \code{ScoreMatrixListControl} function results
#'
#' The resulting object is an extension of a \code{ScoreMatrixList} object, that stores \code{ScoreMatrixControl} objects
#' 
#' @name ScoreMatrixListControl-class
#' @rdname ScoreMatrixListControl-class
#'  
#' @section Constructors:
#' see \code{\link{ScoreMatrixListControl}}
#' 
#' @section Coercion:
#' as(from, "ScoreMatrixListControl"): Creates a \code{ScoreMatrixListControl} object from a list of \code{\link{ScoreMatrixControl}} objects. 
#' 
#' @section Subsetting:
#' In the code snippets below, x is a ScoreMatrixListControl object.
#'  
#' \code{x[[i]]}: Get or set \code{ScoreMatrixControl} object \code{i}, where \code{i} is a numeric or character vector of length 1.
#' 
#' \code{x$name}: value: Get or set element \code{name}, where \code{name} is a name or character vector of length 1.
#' 
#' @slot names A character vector containing names of all stored \code{ScoreMatrixControl} objects.
#' 
#' @seealso \code{\link{ScoreMatrixControl}}, \code{\link{ScoreMatrixList}}
#' @export
setClass("ScoreMatrixListControl",
         slots = c(names = "character"),
         contains = "ScoreMatrixList"
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
#'                  \item{\code{num.annotation}}{vector}

#'                  \item{\code{num.precedence}}{vector}
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
#'                  \item{\code{num.annotation}}{vector}

#'                  \item{\code{num.precedence}}{vector}
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



