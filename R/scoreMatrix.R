#######################################
# S3 functions
#######################################

# returns a matrix
_make.scoreMatrix<-function(target,windows, strand.aware = FALSE){
  
              
            #check if there are common chromsomes
            win.list=as(windows, "RangesList")
            chrs.t=names(target)
            chrs.w=names(win.list)
            chrs  =chrs.t[chrs.t %in% chrs.w]
            if(length(chrs)==0){
              stop("There are no common chromosomes/spaces to do overlap")
            }
            
            
            #get views           
            my.vList=Views( target[names(target) %in% chrs], win.list[names(win.list) %in% chrs] )            
            # get vectors from Views and make a matrix outof it
            my.func<-function(x) t(viewApply( x,as.vector,simplify=TRUE))
            
            # make a matrix from those views            
            mat.list = sapply(my.vList,my.func,simplify=FALSE,USE.NAMES = FALSE);class(mat.list)
            mat = do.call("rbind",mat.list)
			
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



#' get scores overlapping with windows in a scoreMatrix object
#'
#' A scoreMatrix object can be used to draw average profiles or heatmap of read coverage or wig track-like data.
#' \code{windows} can be a predefined region around transcription start sites or other regions of interest that have equal lengths
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
            constraint=GRanges(seqnames=names(r.chr.len),IRanges(start=rep(1,length(r.chr.len)),end=unlist(r.chr.len))  )
            windows=subsetByOverlaps(windows, constraint,type = "within",ignore.strand = TRUE)
            
            mat = _make.scoreMatrix(target, windows, strand.aware=strand.aware)
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



