#######################################
# S3 functions
#######################################


# removes ranges that fell of the rle object
# does not check for the correspondence of the chromosome names - always chech before using this function
constrainRanges = function(target, windows){
	
	checkClass(target, 'SimpleRleList')
	checkClass(windows, 'GRanges')
	
	IRanges::values(windows)$X_rank = 1:length(windows)
	r.chr.len = elementLengths(target)
    constraint = GRanges(seqnames=names(r.chr.len),IRanges(start=rep(1,length(r.chr.len)),end=as.numeric(r.chr.len)))
	# suppressWarnings is done becuause GenomicRanges function give warnings if you don't have the same seqnames in both objects
    win.list.chr = suppressWarnings(IRanges::subsetByOverlaps(windows, constraint,type = "within",ignore.strand = TRUE))
	
	if(length(win.list.chr) == 0)
		stop('All windows fell have coordinates outside chromosome boundaries')
	return(win.list.chr)
}

# given a RleList and a granges object it selcets the intersecting chromosomes and fetches the views
getViews = function(target, windows){

	checkClass(target, 'SimpleRleList')
	checkClass(windows, 'GRanges')

	# orders the granges object so that we can track which view corresponds to which range
	windows = windows[order(as.character(seqnames(windows)), start(windows))]
	win.list=as(windows, "RangesList")
	#check if there are common chromsomes
	chrs  = intersect(names(win.list), names(target))
	if(length(chrs)==0)
		stop("There are no common chromosomes/spaces to do overlap")
		
	#get views 
	# the subsetting needs to be done using a character vector, because otherwise it can take the views from wrong seqnames
	my.vList = Views(target[chrs], win.list[chrs] )
	
	# rownames of each view correspond to the ids of each window
	my.vList = RleViewsList(lapply(chrs, 
								   function(x){
									v = my.vList[[x]]
									names(v) = IRanges::values(windows)$X_rank[as.character(seqnames(windows)) == x]
									v}))
	names(my.vList) = chrs
	return(my.vList)
}

# checkw whether the x object corresponds to the given class
checkClass = function(x, class.name, var.name = deparse(substitute(x))){

	fun.name = match.call(call=sys.call(sys.parent(n=1)))[[1]]
	if(!class(x) == class.name)
		stop(paste(fun.name,': ', var.name, ' is not of class: ', class.name, '\n', sep=''))
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
#' @usage scoreMatrix(target,windows,strand.aware=FALSE,...)
#' @return returns a \code{scoreMatrix} object
#' @seealso \code{\link{scoreMatrixBin}}, \code{\link{modCoverage}}

#' @docType methods
#' @rdname scoreMatrix-methods           
#' @export
setGeneric("scoreMatrix",function(target,windows,strand.aware=FALSE,...) standardGeneric("scoreMatrix") )


# ------------------------------------------------------------------------------------ #
#' @aliases scoreMatrix,RleList,GRanges-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("RleList","GRanges"),
          function(target,windows,strand.aware){
            
            #check if all windows are equal length
            if( length(unique(width(windows))) >1 )
              stop("width of 'windows' are not equal, provide 'windows' with equal widths")
            
			# set a uniq id for the GRanges
			windows = constrainRanges(target, windows)
			
			# fetches the windows
			viewsList = getViews(target, windows)
			mat = do.call(rbind, lapply(viewsList, function(x)t(viewApply(x, as.vector))))
			rownames(mat) = unlist(lapply(viewsList, names), use.names=F)
	
			if(strand.aware == TRUE){
					 #s.ind = as.vector(strand(windows) == '-') # this commented out part will not work, the order of mat is not same as windows
                                         orig.rows=which(as.character(strand(windows))== '-') # this will work
        
                                         mat[rownames(mat) %in% orig.rows,] = mat[rownames(mat) %in% orig.rows, ncol(mat):1]
					#mat[s.ind,] = t(apply(mat[s.ind,],1, rev))
			}
            return(new("scoreMatrix",mat))
})

#' @aliases scoreMatrix,modRleList,GRanges-method
#' @rdname scoreMatrix-methods
setMethod("scoreMatrix",signature("modRleList","GRanges"),
          function(target,windows,strand.aware){
		  
		  mat = scoreMatrix(as(target, 'RleList'), windows, strand.aware)
		  mat=(mat-target@add)/target@multiply
		  mat[mat<0]=NA
		  return(mat)	
})

#' @aliases scoreMatrix,GRanges,GRanges-method
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

#' @param mat a \code{scoreMatrix} object
#' @param fact a \code{factor} of length equal to \code{nrow(mat)}. Unused factor levels are dropped
#' @param ord.vec a \code{vector} of class \code{numeric} of the same length as mat, which is going to be used for ordering of the rows
#' @param shift shift the start coordinate of the x axis (plot starts at -shift)
#' @param mat.cols a vector of colors used for plotting of the heatmap. Default colors range from lightgray to darkblue.
#' @param fact.cols a vector of colors used for plotting of the factor key
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main plot name
#' @param class.names names for each factor class - has to have the same lenght as \code{levels(fact)}
#' @param ... other options to be passed to functions (obsolete at the moment)


#' @usage plotMatrix(mat, fact=NULL, ord.vec=NULL, shift=0, mat.cols=NULL, fact.cols=NULL, xlab='Position', ylab='Region', main='Positional profile', class.names=NULL, ...)
#' @return nothing

#' @docType methods
#' @rdname plotMatrix-methods
#' @export
setGeneric("plotMatrix", function(mat, fact=NULL, ord.vec=NULL, shift=0, mat.cols=NULL, fact.cols=NULL, xlab='Position', ylab='Region', main='Positional profile', class.names=NULL, ...) standardGeneric("plotMatrix") )

#' @aliases plotMatrix,scoreMatrix-method
#' @rdname plotMatrix-methods
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

