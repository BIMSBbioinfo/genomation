# ---------------------------------------------------------------------------#
#' Construct a list of scoreMatrixObjects that can be used for plotting
#'
#' @usage scoreMatrixBin(...), where ... corresponds to 1 or more scoreMatrix objects
#' @return returns a \code{scoreMatrixList} object
#' @export
#' @docType methods
#' @rdname scoreMatrixList-methods        
setGeneric("scoreMatrixList",function(...) standardGeneric("scoreMatrixList"))

# constructor
#' @aliases scoreMatrixList
#' @rdname scoreMatrixList-methods
scoreMatrixList = function(...){
	l = list(...)
	len = length(l)
	if(len == 0L){
		stop('length of the arguments is zero')
	}else{
		if(len == 1L && is.list(l[[1]]))
			l = l[[1]]
		if(!all(unlist(lapply(l, class)) == 'scoreMatrix'))
			stop('Not all provided elements are of class scoreMatrix')
	}
	return(new("scoreMatrixList",l))
}

# --------------------------------- #
# Validator
.valid.scoreMatrixList = function(x){
	errors = character()
	# checks whether all matrices are of class scoreMatrix
	if(!all(unlist(lapply(l, function(x)class(x) == 'scoreMatrix'))))
		errors = paste(errors, 'All elements for scoreMatrixList need to be of class scoreMatrix', sep='\n')

	# checks whether all matrices are numeric
	if(!all(unlist(lapply(l, function(x)all(is.integer(x) | is.numeric(x))))))
		errors = paste(errors, 'Not all matrices are of type integer or numeric', sep='\n')
		
	if(length(errors) == 0) TRUE else errors 
}


# --------------------------------- #
# show Methods
#' @rdname scoreMatrixList-methods
setMethod("show", "scoreMatrixList",
			function(object){
				dims = lapply(object, dim)
				len = length(object)
				widths = apply(do.call(rbind, dims),2, function(x)max(nchar(x)))
				cat('scoreMatrixlist of length:', len, '\n')
				for(i in 1:len){
					s=sprintf(paste('%d%s ','%',widths[1],'d %',widths[2],'d\n', sep=''), 
							i, '. scoreMatrix with dims:', dims[[i]][1], dims[[i]][2])
					cat(s)
				}
			}
)

# --------------------------------- #
# plot functions for score matrix list
#' Plot a scoreMatrixList object, along with the underlying factors
#'
#' @usage heatmapProfile(mat.list, mat.cols)
#' @param'\code{mat.list} a scoreMatrixList object
#' @param'\code{mat.cols} colors to be used for plotting

#' @return returns a \code{scoreMatrixList} object
#' @export
#' @docType methods
#' @rdname heatmapProfile-methods
setGeneric("heatmapProfile", function(mat.list, mat.cols=NULL, ...)standardGeneric("heatmapProfile"))

#' @rdname heatmapProfile-methods
setMethod("heatmapProfile", signature("scoreMatrixList"),
			function(mat.list, mat.cols){
				
				dims = unlist(lapply(mat.list, nrow))
				if(!length(unique(dims)) == 1)
					stop('scoreMatrixList does not contain matrices with the same number of rows')
				
				if(is.null(mat.cols)){
					cat('Using default mat.cols...\n')
					mat.cols = colorRampPalette(c('lightgray','darkblue'), interpolate='spline')(20)
				}
				
				ncols = unlist(lapply(mat.list, ncol))
				nrow = nrow(mat.list[[1]])
				len = length(mat.list)
				layout(matrix(1:len, ncol=len), ncols)

				for(i in 1:len){
					cat('Plotting matrix: ', i,'\r')
					if(i == 1){
						par(mar=c(3,3,3,.5))
					}else{
						par(mar=c(3,.5,3,.5))
					}
					
					image(x=1:ncols[1], y=1:nrow, z=t(mat.list[[i]]), main=names(mat.list)[i], col=mat.cols,yaxt='n', xlab='', ylab='', useRaster=T)
					if(i==1){
						s = round(fivenum(1:nrow))
						axis(2, at=s, las=2, cex.lab=2, labels=s)
					}
				}
				cat('n')
		 }
)

