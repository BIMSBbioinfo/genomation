
# ----------------------------------------------------------------------------------------------- #
#' Gets a weighted coverage vector from a GRanges object with a numeric column
#'
#' This function is similar to \code{coverage} function of IRanges and GenomicRanges, however it can deal with floating point numbers to some degree.
#' The resulting RleList has values from the denoted numeric column, it is not simply a coverage vector. In order to deal with numbers with decimals some heuristics  are applied.
#' Warning: Users should NOT use this unless they understand what it exactly does.The default arguments will just return a coverage vector weighted by the denoted column of the GRanges object
#'
#' @param x GRanges object with elementData
#' @param col.name name of the column that has a numeric value to be used as a value in RleList object. It should be a column from the \code{elementMetadata} section of the object.
#' @param multiply a value (default:1000) that will be multiplied by the value denoted by \code{col.name} in GRanges object. This is useful to retain some of the values with decimals.
#' @param add This value (default:1) be added to the value denoted by \code{col.name} in GRanges object. This addition will procede the multiplication. This way the bases that have no defined value will have 0 value in the Rle vector, and the bases with a predefined value (even though it is 0 will be distinguishable in subsequent operations)
#'
#' @usage modCoverage(x,col.name,multiply=1000,add=1)
#' @return returns a \code{modRleList} object
#'
#' @note \code{multiply} argument helps users retain some of the values with decimal points. Users have to keep in mind the fact that \code{add} and/or \code{multiply} arguments are in effect in the resulting RleList object and do subsequent calculations keep that fact in mind.
#'
#'
#' @seealso \code{\link{coverage}},\code{\link{modRleList}}
#' @export
#' @docType methods
#' @rdname modCoverage-methods
setGeneric("modCoverage", function(x,col.name,multiply=1000,add=1) standardGeneric("modCoverage") )

#' @aliases modCoverage,GRanges,character-method
#' @rdname modCoverage-methods
setMethod("modCoverage",
			signature(x = "GRanges", col.name = "character"),
			function( x,col.name,multiply,add ){
				if(  !is.numeric(elementMetadata(x)[col.name][,1]) )
					stop("The column you provided ", col.name," is non-numeric! Specify a numeric column")

				new("modRleList",
					coverage(x, weight = (multiply*elementMetadata(x)[col.name][,1])+add ),
					multiply=multiply,
					add=add)
})
