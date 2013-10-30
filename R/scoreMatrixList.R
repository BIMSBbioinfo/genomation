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

# ---------------------------------------------------------------------------- #
# ScoreMatrixList constructor
#' Construct a list of scoreMatrixObjects that can be used for plotting
#'
#' @param l: corresponds to can be a list of \code{scoreMatrix} objects, that are coerced to the \code{ScoreMatrixList}, or a list of \code{RleList} objects that are used to construct the \code{scoreMatrixList}
#' @param granges: a \code{GenomicRanges} containing viewpoints for the scoreMatrix or ScoreMatrixList functions
#' @param bin: an integer telling the number of bins to bin the score matrix
 
#' @usage ScoreMatrixList(l, granges, bin)
#' @return returns a \code{ScoreMatrixList} object
#' @export
#' @docType methods
#' @rdname ScoreMatrixList-methods
ScoreMatrixList = function(l, granges=NULL, bin=NULL, ...){

	len = length(l)
	if(len == 0L)
		stop('list argument is empty')
	
	# checks whether the list argument contains only scoreMatrix objects
	if(all(unlist(lapply(l, class)) == 'scoreMatrix'))
		return(new("ScoreMatrixList",l))
	
	# Given a list of RleList objects and a granges object, returns the scoreMatrix list Object
	if(all(unlist(lapply(l, class)) %in% c('SimpleRleList', 'RleList'))){
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
	
	return(new("ScoreMatrixList",sml))
}

# ---------------------------------------------------------------------------- #
# Validator
.valid.ScoreMatrixList = function(x){
	errors = character()
	# checks whether all matrices are of class scoreMatrix
	if(!all(unlist(lapply(l, function(x)class(x) == 'scoreMatrix'))))
		errors = paste(errors, 
                   'All elements for ScoreMatrixList 
                    need to be of class scoreMatrix', 
                    sep='\n')

	# checks whether all matrices are numeric
	if(!all(unlist(lapply(l, function(x)all(is.integer(x) | is.numeric(x))))))
		errors = paste(errors, '
                   Not all matrices are of type integer or numeric', 
                   sep='\n')
		
	if(length(errors) == 0) TRUE else errors 
}


# ---------------------------------------------------------------------------- #
# show Methods
#' @rdname show-methods
setMethod("show", "ScoreMatrixList",
			function(object){
				dims = lapply(object, dim)
				len = length(object)
				widths = apply(do.call(rbind, dims),2, function(x)max(nchar(x)))
				message('scoreMatrixlist of length:', len, '\n')
				for(i in 1:len){
					s=sprintf(paste('%d%s ','%',widths[1],'d %',widths[2],'d\n', sep=''), 
							i, '. scoreMatrix with dims:', dims[[i]][1], dims[[i]][2])
					message(s)
				}
			}
)

# ---------------------------------------------------------------------------- #
#' plot functions for score matrix list
#' Plot a ScoreMatrixList object as a panel of heatmaps
#'
#' @param mat.list a \code{ScoreMatrixList} object
#' @param mat.cols colors to be used for plotting
#' @param xmarks an integer number to lable the thick marks on the x axis of each heatmap. By default it takes the values of -ncol/2, 0, ncol/2
#' @param ymarks a vector of that will lable the thick marks on the y axis
#' @param y.at a numeric vector that will specify the positions of the thick marks on the y axis
#' @param xcex, ycex an integer number which controls the character expansion on x and y axis
#' @param cex.main an integer number which controls the character expansion of the plot label
#' @param mar a vector of length 5 which controls the size of the margins. The order is the following: below, left, up, right, spacing between consecutive plots
#' @param use.names whether to use the names of the ScoreMatrixList object to label each plot
#' @param main whether to use the names of the ScoreMatrixList object to label each plot
#' @param xlab, ylab name to be used for the x/y axis 
#' @param ... other options (obselete for now)

#' @usage heatmapProfile(mat.list, mat.cols=NULL, ...)
#' @docType methods
#' @rdname ScoreMatrixList-methods

#' @example l = lapply(seq(20, 40,5), function(x)as(matrix(rpois(1000, x), ncol=25), 'scoreMatrix'))
#' @example l = ScoreMatrixList(l)
#' @example names(l) = letters[1:5]
#' @example heatmapProfile(l)
#' @export
setGeneric("heatmapProfile", 
           function(mat.list, mat.cols=NULL, xmarks=NULL, 
                    ymarks=NULL, y.at=NULL, xcex=1.5, 
                    ycex=1.5, cex.main=3, mar=NULL, 
                    use.names=T, xlab=NULL, ylab=NULL, ...)
                      standardGeneric("heatmapProfile"))

#' @aliases heatmapProfile,ScoreMatrixList-method
#' @rdname heatmapProfile-methods
setMethod("heatmapProfile", signature(mat.list="ScoreMatrixList"),
			function(mat.list, mat.cols, xmarks, 
               ymarks, y.at, xcex, ycex, 
               cex.main, mar, use.names, xlab, ylab, ...){
				
				dims = unlist(lapply(mat.list, nrow))
				if(!length(unique(dims)) == 1)
					stop('ScoreMatrixList does not contain matrices with the same number of rows')
								
				
				# default matrix colors
				if(is.null(mat.cols))
					mat.cols = colorRampPalette(c('lightgray','darkblue'), 
                                      interpolate='spline')(20)

				
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
					layout(matrix(c(2:(len+1), rep(1, len)), 
                        ncol=len, nrow=2, byrow=T), height=c(20,1))
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
				
				# sets the thick marks and labels on the y axis
				if(is.null(ymarks)){
					ymarks = round(fivenum(1:nrow))
				}
				if(is.null(y.at)){
					y.at = round(fivenum(1:nrow))
				}
				if(length(ymarks) != length(y.at))
					stop('ymarks and y.at do not have the same length')
				if(!is.null(y.at)){
					if(!is.numeric(y.at)){
						stop('y.at need to be a numeric variable')
					}
					if(any((y.at < 1) | (y.at > nrow))){
						stop('y.at values are outside of the matrix dimension')
					}
				}
				
				# plots the heatmaps
				for(i in 1:len){
					message('Plotting matrix: ', i,'\r')
					# sets the margins for each plot
					if(i == 1){	
						par(mar=c(mar[c(1,2,3)], mar[5]/2))

					}else if(i == len){
						par(mar=c(mar[1], mar[5]/2, mar[3], mar[4]))
					}else{
						par(mar=c(mar[1], mar[5]/2, mar[3], mar[5]/2))
					}
					
					image(x=1:ncols[1], y=1:nrow, z=t(mat.list[[i]]), 
                main=main[i], cex.main=cex.main, col=mat.cols,
                yaxt='n', xaxt='n', xlab='', ylab='', useRaster=T, ...)
					if(i==1){
						axis(2, at=y.at, las=2, cex.lab=2, labels=ymarks, cex.axis=ycex)
					}
					
					axis(1, at=xpos, labels=xmarks, cex.axis=xcex)
				}
				message('\nPlotting done\n')
				
		 }
)


# ---------------------------------------------------------------------------- #
#' Scales each scoreMatrix in the ScoreMatrixList object

#' @param sml a \code{ScoreMatrixList} object
#' @param columns a \code{columns} whether to scale the matrix by columns. Set by default to FALSE.
#' @param rows  a \code{rows} Whether to scale the matrix by rows. Set by default to TRUE
#' @param scalefun a function object that takes as input a matrix and returns a matrix. By default  the argument is set to the R scale function with center=TRUE and scale=TRUE

#' @usage scaleScoreMatrixList(mat, columns=FALSE, rows=TRUE, ...)
#' @return \code(ScoreMatrixList) object

#' @docType methods
#' @rdname ScoreMatrixList-methods
#' @export
setGeneric("scaleScoreMatrixList", 
           function(sml, 
                    columns=FALSE, rows=TRUE, 
                    scalefun=function(x)scale(x), 
                    ...) 
             standardGeneric("scaleScoreMatrixList") )

setMethod("scaleScoreMatrixList", signature("ScoreMatrixList"),
          function(sml, columns, rows, scalefun, ...){
            
            sml = lapply(sml, function(x)
                          scaleScoreMatrix(x, 
                                           columns=colums, 
                                           rows=rows, 
                                           scalefun=scalefun))
            sml = ScoreMatrixList(sml)
            return (sml)
          }
)

# ---------------------------------------------------------------------------- #
#' Returns a union of rows for each matrix in a ScoreMatrixList object. 
#' This is done using the rownames of each element in the list.

#' @param sml a \code{ScoreMatrixList} object


#' @usage unionScoreMatrixList(sml, columns=FALSE, rows=TRUE, ...)
#' @return \code(ScoreMatrixList) object

#' @docType methods
#' @rdname ScoreMatrixList-methods
#' @export
setGeneric("unionScoreMatrixList", 
           function(sml)
             standardGeneric("unionScoreMatrixList") )

setMethod("unionScoreMatrixList", signature("ScoreMatrixList"),
          function(sml){
            
            rnames = Reduce('union' ,lapply(sml, rownames))
            sml = ScoreMatrixList(
                    lapply(sml, function(x)x[rownames(x) %in% rnames,]))
            return (sml)
          }
)