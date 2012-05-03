#######################
#         S3          #
#######################

.plotXYlab = function(lab, pos){
	
	if(!pos %in% c('x','y'))
		stop('pos can only be x or y')
	
	if(pos == 'x'){
		.plotLab(lab, srt=0)
	}

	if(pos == 'y'){
		.plotLab(lab, srt=90)
	}
}

.plotLab = function(lab, srt, xlim, ylim, x, y){
	plot.new()
	plot.window(xlim=c(1,5), ylim=c(1,5))
	text(3, 3, labels=lab, cex=2, srt=srt)
}

# ------------------------------------------------------------------------------------ #
# scoreMatrixList constructor
#' Construct a list of scoreMatrixObjects that can be used for plotting
#'
#' @param l: corresponds to can be a list of \code{scoreMatrix} objects, that are coerced to the \code{scoreMatrixList}, or a list of \code{RleList} objects that are used to construct the \code{scoreMatrixList}
#' @param granges: a \code{GenomicRanges} containing viewpoints for the scoreMatrix or scoreMatrixList functions
#' @param bin: an integer telling the number of bins to bin the score matrix
 
#' @usage scoreMatrixList(l, granges, bin)
#' @return returns a \code{scoreMatrixList} object
#' @export
#' @docType methods
#' @rdname scoreMatrixList-methods
scoreMatrixList = function(l, granges=NULL, bin=NULL, ...){

	len = length(l)
	if(len == 0L){
		stop('list argument is empty')
	
	# checks whether the list argument contains only scoreMatrix objects
	if(all(unlist(lapply(l, class)) == 'scoreMatrix'))
		return(new("scoreMatrixList",l))
	
	# Given a list of RleList objects and a granges object, returns the scoreMatrix list Object
	if(all(unlist(lapply(l, class)) %in% c('SimpleRleList', 'RleList')){
		if(is.null(granges))
			stop("If the list contains RleLists granges must be defined")
		
		
		if(is.null(bin) && all(width(granges)) == unique(width(granges))){
			sml = lapply(l, function(x)scoreMatrix(x, granges))
		
		} else{
			if(is.null(bin))
				bin = 10
			sml = lapply(l, function(x)scoreMatrixBin(x, granges, bin=bin))
		}
	
	}else{
		stop("List does not contain the proper classes")
	}
	
	return(new("scoreMatrixList",sml))
}

# ------------------------------------------------------------------------------------ #
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


# ------------------------------------------------------------------------------------ #
# show Methods
#' @rdname show-methods
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

# ------------------------------------------------------------------------------------ #
#' plot functions for score matrix list
#' Plot a scoreMatrixList object as a panel of heatmaps
#'
#' @param mat.list a \code{scoreMatrixList} object
#' @param mat.cols colors to be used for plotting
#' @param xmarks an integer number to lable the thick marks on the x axis of each heatmap. By default it takes the values of -ncol/2, 0, ncol/2
#' @param xcex, ycex an integer number which controls the character expansion on x and y axis
#' @param cex.main an integer number which controls the character expansion of the plot label
#' @param mar a vector of length 5 which controls the size of the margins. The order is the following: below, left, up, right, spacing between consecutive plots
#' @param use.names whether to use the names of the scoreMatrixList object to label each plot
#' @param main whether to use the names of the scoreMatrixList object to label each plot
#' @param xlab, ylab name to be used for the x/y axis 
#' @param ... other options (obselete for now)

#' @usage heatmapProfile(mat.list, mat.cols=NULL, ...)
#' @docType methods
#' @rdname scoreMatrixList-methods

#' @example l = lapply(seq(20, 40,5), function(x)as(matrix(rpois(1000, x), ncol=25), 'scoreMatrix'))
#' @example l = scoreMatrixList(l)
#' @example names(l) = letters[1:5]
#' @example heatmapProfile(l)
#' @export
setGeneric("heatmapProfile", function(mat.list, mat.cols=NULL, xmarks=NULL, xcex=1.5, ycex=1.5, cex.main=3, mar=NULL, use.names=T, xlab=NULL, ylab=NULL, ...)standardGeneric("heatmapProfile"))

#' @aliases heatmapProfile,scoreMatrixList-method
#' @rdname heatmapProfile-methods
setMethod("heatmapProfile", signature(mat.list="scoreMatrixList"),
			function(mat.list, mat.cols, xmarks, xcex, ycex, cex.main, mar, use.names, xlab, ylab, ...){
				
				dims = unlist(lapply(mat.list, nrow))
				if(!length(unique(dims)) == 1)
					stop('scoreMatrixList does not contain matrices with the same number of rows')
								
				
				# default matrix colors
				if(is.null(mat.cols))
					mat.cols = colorRampPalette(c('lightgray','darkblue'), interpolate='spline')(20)

				
				# checks the margin parameter
				if(is.null(mar)){
					mar = rep(3, 5)
				}else if(! length(mar) == 5 | !all(is.numeric(mar))){
					stop('mar is not of length 5')
				}
				
				if(use.names)
					main = names(mat.list)
				
				# gets the dimension of the matrices
				ncols = unlist(lapply(mat.list, ncol))
				nrow = nrow(mat.list[[1]])
				len = length(mat.list)
				
				# sets the layout
				if(is.null(xlab) & is.null(ylab)){
					layout(matrix(1:len, ncol=len))
				
				}else if(! is.null(xlab) & is.null(ylab)){
					layout(matrix(c(2:(len+1), rep(1, len)), ncol=len, nrow=2, byrow=T), height=c(20,1))
					.plotXYlab(xlab, 'x')
				
				}else if(is.null(xlab) & ! is.null(ylab)){
					layout(matrix(1:(len+1), ncol=len+1), width=c(5,rep(20, len+1)))
					.plotXYlab(ylab, 'y')
				
				}else if(!is.null(xlab) & !is.null(ylab)){
					par(mar=c(0,1,0,1), oma=rep(0,4))
					layout(matrix(c(1,3:(len+2), 0,rep(2, len)), ncol=len+1, byrow=T), 
						   width=c(4,rep(20, len+1)), 
						   height=c(20,1))
					.plotXYlab(ylab, 'y')
					.plotXYlab(xlab, 'x')
				}
				
				# gets the tick marks labels and positions
				if(is.null(xmarks))
					xmarks = c(-ncols[1]/2, 0, ncols[1]/2)
				xpos = seq(1, ncols[1], length.out=length(xmarks))
				
				for(i in 1:len){
					cat('Plotting matrix: ', i,'\r')
					# sets the margins for each plot
					if(i == 1){
						par(mar=c(mar[c(1,2,3)], mar[5]/2))
					}else if(i == len){
						par(mar=c(mar[1], mar[5]/2, mar[3], mar[4]))
					}else{
						par(mar=c(mar[1], mar[5]/2, mar[3], mar[5]/2))
					}
					
					image(x=1:ncols[1], y=1:nrow, z=t(mat.list[[i]]), main=main[i], cex.main=cex.main, col=mat.cols, yaxt='n', xaxt='n', xlab='', ylab='', useRaster=T, ...)
					if(i==1){
						s = round(fivenum(1:nrow))
						axis(2, at=s, las=2, cex.lab=2, labels=s, cex.axis=ycex)
					}
					
					axis(1, at=xpos, labels=xmarks, cex.axis=xcex)
				}
				cat('\nPlotting done\n')
				
		 }
)

