#######################################
# S3 functions
#######################################

# returns a matrix
.make.scoreMatrix<-function(target,windows, strand.aware = FALSE){
  
             
            
            #get views           
            my.vList = Views(target[chrs], win.list[chrs] )            
            # get vectors from Views and make a matrix outof it
            my.func = function(x) t(viewApply( x, as.vector,simplify=TRUE))
            
            # make a matrix from those views            
            mat.list = sapply(my.vList,my.func,simplify=FALSE,USE.NAMES = FALSE)
            mat = do.call("rbind", mat.list)
			rownames(mat) = 
			
			# if the order is strand aware it reverses the profiles on the negative strand
			if(strand.aware == TRUE){
				s.ind = as.vector(strand(windows) == '-')
				mat[s.ind] = rev(mat[s.ind])
			}
                    
            # if the target is modRleList do appropriate calculations to get the score and put NAs in cells that have no value
            if(is(target,"modRleList")){
              
              #remove values full of NA values
              #mat=mat[rowSums(mat)>0,]
              
              mat=(mat-target@add)/target@multiply
              mat[mat<0]=NA
              
            }
      return(mat)
}


#######################################
# S4 functions
#######################################


# ------------------------------------------------------------------------------------ #
#' get scores overlapping with windows in a scoreMatrix object
#'
#' A scoreMatrix object can be used to draw average profiles or heatmap of read coverage or wig track-like data.
#' \code{windows} can be a predefined region around transcription start sites or other regions of interest that have equal lengths
#' The function removes all window that fall off the Rle object - have the start coordinate < 1 or end coordinate > length(Rle)
#' The function takes the intersection of names in the Rle and GRanges objects
#'
#' @param target a \code{RleList} or a \code{modRleList} or \code{GRanges} object to be overlapped with ranges in \code{windows}
#' @param windows a \code{GRanges} object that contains the windows of interest. It could be promoters, CpG islands, exons, introns. However the sizes of windows have to be equal.
#' @param strand.aware If TRUE (default: FALSE), the strands of the windows will be taken into account in the resulting \code{scoreMatrix}. If the strand of a window is -, the values of the bins for that window will be reversed
#' @param ... parameters to be passed to \code{modCoverage} function. Only needed when target is \code{GRanges}.
#'
#' @usage scoreMatrix(target,windows,target,windows,strand.aware=FALSE,...)
#' @return returns a \code{scoreMatrix} object
#' @seealso \code{\link{scoreMatrixBin}}, \code{\link{modCoverage}}

#' @export
#' @docType methods
#' @rdname scoreMatrix-methods           
setGeneric("scoreMatrix",function(target,windows,strand.aware=FALSE,...) standardGeneric("scoreMatrix") )


# ------------------------------------------------------------------------------------ #
#' @aliases scoreMatrix,GRanges,RleList-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("RleList","GRanges"),
          function(target,windows,strand.aware){
            
            #check if all windows are equal length
            if( length(unique(width(windows))) >1 ){
              stop("width of 'windows' are not equal, provide 'windows' with equal widths")
            }
			
            #check if windows lengths exceeds the length of feature based chromosomes
            r.chr.len=lapply(target,length)
            constraint=GRanges(seqnames=names(r.chr.len),IRanges(start=rep(1,length(r.chr.len)),end=unlist(r.chr.len)))
			values(windows)$rank = 1:length(windows)
            windows=subsetByOverlaps(windows, constraint,type = "within",ignore.strand = TRUE)
            
            mat = .make.scoreMatrix(target, windows, strand.aware=strand.aware)
            return( new("scoreMatrix",mat) )
})

#' @aliases scoreMatrix,GRanges,GRanges,ANY-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("GRanges","GRanges"),
          function(target,windows,strand.aware,...){
            
            #make coverage vector (modRleList) from target
            target.rle=modCoverage(target,...)
            
            # call scoreMatrix function
            scoreMatrix(target.rle,windows,strand.aware)
})


# ------------------------------------------------------------------------------------ #
#' visual representation of scoreMatrix using a heatmap 
#' The rows can be reordered using one factor and one numeric vector

#' @param mat a \code(scoreMatrix) object
#' @param fact a \code(factor) of length equal to \code(nrow(mat)). Unused factor levels are dropped
#' @param ord.vec a \code(vector) of class \code(numeric) of the same length as mat, which is going to be used for ordering of the rows
#' @param shift shift the start coordinate of the x axis (plot starts at -shift)
#' @param mat.cols a vector of colors used for plotting of the heatmap. Default colors range from lightgray to darkblue.
#' @param fact.cols a vector of colors used for plotting of the factor key
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main plot name
#' @param class.names names for each factor class - has to have the same lenght as \code(levels(fact))


#' @usage plotMatrix(mat, fact, ord.vec, shift, mat.cols, ord.vec, shift, mat.cols, fact.cols, xlab, ylab, main, ...)
#' @return nothing

#' @seealso
#' @docType methods
#' @rdname scoreMatrix-methods
setGeneric("plotMatrix", function(mat, fact=NULL, ord.vec=NULL, shift=0, mat.cols=NULL, fact.cols=NULL, xlab='Position', ylab='Region', main='Positional profile', class.names=NULL, ...) standardGeneric("plotMatrix") )

#' @rdname scoreMatrix-methods
setMethod("plotMatrix", signature("scoreMatrix"),
		  function(mat, fact, ord.vec, shift, mat.cols, fact.cols, xlab, ylab, main, class.names, ...){
			
			# -------------------------- #
			# parameter checking
			if(!is.null(fact) & length(fact) != nrow(mat))
				stop('Given factor does not have the same length as the matrix')
			if(!is.null(fact) & !is.factor(fact))
				stop('fact need to be an object of class factor')
			if(!is.null(ord.vec) & length(fact) != nrow(mat))
				stop('Given ordering vector does not have the same length as the matrix')
			if(!is.numeric(shift) | length(shift) > 1)
				stop('shift needs to be a numeric vector of length 1')
			# -------------------------- #
			# default values
			if(is.null(fact))
				fact=as.factor(rep(1, nrow(mat)))
			if(is.null(ord.vec))
				ord.vec = 1:nrow(mat)	
			if(is.null(class.names))
				class.names = levels(fact)
			
			if(is.null(mat.cols)){
				cat('Using default mat.cols...\n')
				mat.cols = colorRampPalette(c('lightgray','darkblue'), interpolate='spline')(20)
			}
			if(is.null(fact.cols)){
				cat('Using default fact.cols...\n')
				fact.cols = adjustcolor(rainbow(length(levels(fact))), offset=c(0.5,0.5,0.5, 0), transform=diag(c(.7, .7, .7, 0.6)))
				 
			}
			
			# drops unused levels from the factor
			fact = fact[1:length(fact), drop=TRUE]
			
			# -------------------------- #
			# plots the matrix
			mat = mat[order(as.numeric(fact), ord.vec),]
			# par(fig=c(0,.95,0,1), mar=c(5,5,3,.5))
			layout(matrix(c(1,2), ncol=2), widths=c(10,1))
			par(mar=c(5,5,3,.5), oma=c(0,0,0,0))
			AddSep = function(x, rowsep, col, sepwidth=c(0.05,0.5)){
				for(rsep in rowsep){
					rect(xleft =0, ybottom= (rsep), xright=ncol(x)+1,  ytop = (rsep+1) - sepwidth[2], lty=1, lwd=1, col=col, border=col)
				}
			}
			# plots the main matrix
			image(x=1:ncol(mat) - shift, y=1:nrow(mat), z=t(as.matrix(mat)), col=mat.cols, , oma=c(0,0,0,0), useRaster=T, xlab=xlab, ylab=ylab, main=main)
			classnum = table(fact)
			rowsep = cumsum(classnum)
            AddSep(mat, rowsep[-length(rowsep)], "black")	
			
			# plots the class designation
			# par(fig=c(.95,1,0,1), new=TRUE, mar=c(5,.5,3,1))
			par(mar=c(5,.5,3, max(max(nchar(class.names)/2, 5))))
			image(x = 1:20, y = 1:nrow(mat), z=t(matrix(as.numeric(fact), nrow=length(fact), ncol=20)), col = fact.cols, xaxt='n', yaxt='n', ylab='', xlab='')
			at = classnum/2
			at[-1] = at[-1] + at[-length(at)]
			at = cumsum(at)
			axis(side=4, at=at, labels=class.names, tick = F, las=2)
	  }
)

