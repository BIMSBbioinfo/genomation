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
#' ScoreMatrixList constructor
#' 
#' Construct a list of scoreMatrixObjects that can be used for plotting
#'
#' @param l can be a list of \code{scoreMatrix} objects, that are coerced to the \code{ScoreMatrixList}, a list of \code{RleList} objects, or a character vector specifying the locations of mulitple bam files that are used to construct the \code{scoreMatrixList}. If l is either a RleList object or a character vector of files, it is obligatory to give a granges argument.
#' @param granges a \code{GenomicRanges} containing viewpoints for the scoreMatrix or ScoreMatrixList functions
#' @param bin.num an integer telling the number of bins to bin the score matrix
#' @param bin.op an name of the function that will be used for smoothing windows of ranges
#' @param strand.aware a boolean telling the function whether to reverse the coverage of ranges that come from - strand (e.g. when plotting enrichment around transcription start sites)
#' @param type if l is a character vector of file paths, then type designates the type of the corresponding files (bam or bigWig)
#' @param ... other arguments that can be passed to the function
 
#' @return returns a \code{ScoreMatrixList} object
#' @export
#' @docType methods
#' @rdname ScoreMatrixList-methods
ScoreMatrixList = function(target, windows=NULL, bin.num=NULL, 
                           bin.op='mean', strand.aware=FALSE, ...){

	len = length(target)
	if(len == 0L)
		stop('target argument is empty')
  
	# ----------------------------------------------------------------- #
	# checks whether the list argument contains only scoreMatrix objects
	if(all(unlist(lapply(target, class)) == 'scoreMatrix'))
		return(new("ScoreMatrixList",l))

  
	# ----------------------------------------------------------------- #
	if(is.null(windows))
	  stop("windows of class GRanges must be defined")
  
	# Given a list of RleList objects and a granges object, returns the scoreMatrix list Object
	if(!all(unlist(lapply(target, class)) %in% c('SimpleRleList', 'RleList','GRanges')) &
	   !all(file.exists(target)))
      stop('target should be one of the following: 
           an RleList, a list of files, a list of GRanges')
	
  if(all(file.exists(target)) & is.null(type))
      stop('When providing a file, it is necessary to give the type of the file')
  
  sml = list()
  for(i in 1:length(target)){
    
    message(paste('reading file:',basename(target[i])))
    if(is.null(bin.num) && all(width(windows) == unique(width(windows)))){
      sml[[i]] = ScoreMatrix(target[[i]], windows, strand.aware=strand.aware 
                             , ...)
      
    } else{
      if(is.null(bin.num))
        bin.num = 10
      sml[[i]] = ScoreMatrixBin(target[[i]], windows, 
                                bin.num=bin.num, bin.op=bin.op, 
                                strand.aware=strand.aware, ...)
    }  
  }
	
  if(class(target) %in% c('SimpleRleList', 'RleList','GenomicRanges'))
    names(sml) = names(target)
  if(all(is.character(target)))
    names(sml) = basename(target)
  
  
	return(new("ScoreMatrixList",sml))
}

# ---------------------------------------------------------------------------- #
# Validator
.valid.ScoreMatrixList = function(l){
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
					s=sprintf(paste('%d%s ','%',widths[1],'d %',widths[2],'d', sep=''), 
							i, '. scoreMatrix with dims:', dims[[i]][1], dims[[i]][2])
					message(s)
				}
			}
)


# ---------------------------------------------------------------------------- #
#' Scale the ScoreMatrixList
#' 
#' Scales each scoreMatrix in the ScoreMatrixList object, by rows and/or columns
#' 
#' @param sml a \code{ScoreMatrixList} object
#' @param columns a \code{columns} whether to scale the matrix by columns. Set by default to FALSE
#' @param rows  a \code{rows} Whether to scale the matrix by rows. Set by default to TRUE
#' @param scalefun a function object that takes as input a matrix and returns a matrix.
#' @param ... other argments that be passed to the function
#'  By default  the argument is set to the R scale function with center=TRUE and scale=TRUE
#'
#' @usage scaleScoreMatrixList(sml, columns, rows, scalefun, ...)
#' @return \code{ScoreMatrixList} object
#'
#' @docType methods
#' @rdname scaleScoreMatrixList
#' @export
setGeneric("scaleScoreMatrixList", 
           function(sml, 
                    columns=FALSE, rows=TRUE, 
                    scalefun=NULL, 
                    ...) 
             standardGeneric("scaleScoreMatrixList") )

#' @aliases scaleScoreMatrixList
#' @rdname scaleScoreMatrixList
setMethod("scaleScoreMatrixList", signature("ScoreMatrixList"),
          function(sml, columns, rows, scalefun, ...){
            
            
            sml = lapply(sml, function(x)
              scaleScoreMatrix(x, 
                               columns=columns, 
                               rows=rows, 
                               scalefun=scalefun))
            sml = as(sml,'ScoreMatrixList')
            return (sml)
          }
)

# ---------------------------------------------------------------------------- #
#' Get common rows from all matrices in a ScoreMatrixList object
#' 
#' Returns a intersection of rows for each matrix in a ScoreMatrixList object. 
#' This is done using the rownames of each element in the list.
#'
#' @param sml a \code{ScoreMatrixList} object
#' @param reorder if TRUE \code{ScoreMatrix} objects in the list are sorted
#'                based on their common row ids.
#'
#' @return \code{ScoreMatrixList} object
#'
#' @docType methods
#' @rdname intersectScoreMatrixList-methods
#' @export
setGeneric("intersectScoreMatrixList", 
           function(sml,reorder=FALSE)
             standardGeneric("unionScoreMatrixList") )

#' @aliases intersectScoreMatrixList,ScoreMatrixList-method
#' @rdname intersectMatrixList-methods
setMethod("intersectScoreMatrixList", signature("ScoreMatrixList"),
          function(sml,reorder){
            
            rnames = Reduce('intersect' ,lapply(sml, rownames))
            if(reorder){
              sml = as(lapply(sml, function(x){ 
                                  x=x[rownames(x) %in% rnames,]
                                  x[order(rownames(x)),]
                                  }), 
                     'ScoreMatrixList')
            }else{
              sml = as(lapply(sml, function(x)x[rownames(x) %in% rnames,]), 
                       'ScoreMatrixList')              
            }
            return (sml)
          }
)

# ---------------------------------------------------------------------------- #
#' Reorder all elements of a ScoreMatrixList to a given ordering vector

#' @param sml \code{ScoreMatrixList} object
#' @param ord.vec an integer vector
#' @usage order(sml, ord.vec)
#' @return \code{ScoreMatrixList} object

#' @docType methods
#' @export
setMethod("order", signature("ScoreMatrixList"),
          function(sml, ord.vec){
            
            ord.vec = as.integer(ord.vec)
            sml = lapply(sml, function(x)x[ord.vec,])
            return (as(sml,'ScoreMatrixList'))
          }
)


# ---------------------------------------------------------------------------- #
#' @aliases binMatrix,ScoreMatrixList-method

#' @rdname binMatrix-methods
setMethod("binMatrix", signature("ScoreMatrixList"),
          function(x, bin.num=NULL, fun='mean', ...){
            
            if(is.null(bin.num))
              return(x)
             
            return(new("ScoreMatrix", 
                       lapply(x, function(y)binMatrix(y, bin.num=bin.num, fun=fun))))
          }
)
