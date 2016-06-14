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
#' Make ScoreMatrixList from multiple targets
#' 
#' The function constructs a list of \code{ScoreMatrix} objects in the form
#' of \code{ScoreMatrixList} object. This object can be visualized using 
#' \code{multiHeatMatrix}, \code{heatMeta} or \code{plotMeta}
#'
#' @param targets can be a list of \code{scoreMatrix} objects, that are coerced 
#'        to the \code{ScoreMatrixList}, a list of \code{RleList} objects, or a 
#'        character vector specifying the locations of mulitple bam files  or
#'        bigWig files that 
#'        are used to construct the \code{scoreMatrixList}. If it is either a 
#'        RleList object or a character vector of files, it is obligatory to 
#'        give a windows argument.
#' @param windows \code{GenomicRanges} containing viewpoints for the scoreMatrix 
#'        or ScoreMatrixList functions
#' @param bin.num an integer telling the number of bins to bin the score matrix
#' @param bin.op an name of the function that will be used for smoothing windows of ranges
#' @param strand.aware a boolean telling the function whether to reverse the 
#'        coverage of ranges that come from - strand (e.g. when plotting 
#'        enrichment around transcription start sites)
#' @param weight.col if the object is \code{GRanges} object a numeric column
#'                 in meta data part can be used as weights. This is particularly
#'                useful when genomic regions have scores other than their
#'                coverage values, such as percent methylation, conservation
#'                scores, GC content, etc. 
#' @param is.noCovNA (Default:FALSE)
#'                  if TRUE,and if 'targets' is a GRanges object with 'weight.col'
#'                   provided, the bases that are uncovered will be preserved as
#'                   NA in the returned object. This useful for situations where
#'                   you can not have coverage all over the genome, such as CpG
#'                    methylation values.
#'                    
#' @param type if \code{targets} is a character vector of file paths, then type 
#'        designates the type of the corresponding files (bam or bigWig)
#' @param rpm boolean telling whether to normalize the coverage to per milion reads. 
#'            FALSE by default. See \code{library.size}.
#' @param unique boolean which tells the function to remove duplicated reads 
#'                       based on chr, start, end and strand
#' @param extend numeric which tells the function to extend the features
#'               ( i.e aligned reads) to total
#'               length ofwidth+extend
#' @param param ScanBamParam object  
#' @param library.size a numeric vector of the same length as \code{targets} 
#'                     indicating total number of mapped reads in BAM files (\code{targets}).
#'                     If is not given (default: NULL) then library sizes for every target
#'                     is calculated using a Samtools idxstats like function:
#'                     sum(idxStats(target)$mapped).
#'                     \code{rpm} argument has to be set to TRUE.
#' @param cores the number of cores to use (default: 1)
#'
#' @return returns a \code{ScoreMatrixList} object
#' 
#' @examples
#' 
#' # visualize the distribution of cage clusters and cpg islands around promoters
#' library(GenomicRanges)
#' data(cage)
#' data(cpgi)
#' data(promoters)
#'  
#' cage$tpm = NULL
#' targets = GRangesList(cage=cage, cpgi=cpgi)
#' sml = ScoreMatrixList(targets, promoters, bin.num=10, strand.aware=TRUE)
#' sml
#' \donttest{
#' multiHeatMatrix(sml)
#' }
#' @export
#' @docType methods
#' @rdname ScoreMatrixList-methods
ScoreMatrixList = function(targets, windows=NULL, bin.num=NULL, 
                           bin.op='mean', strand.aware=FALSE, weight.col=NULL, 
                           is.noCovNA=FALSE, type='', rpm=FALSE, unique=FALSE, 
                           extend=0, param=NULL, library.size=NULL, cores=1){
  len = length(targets)
  if(len == 0L)
    stop('targets argument is empty')
  
  # this checks whether we can work with the corresponding targets object class set
  list.ind = grepl('list', class(targets)) | grepl('List', class(targets))
  if(len > 1L & !list.ind){
    if(all(is.character(targets)) && is.null(type) && all(file.exists(targets)))
      stop('targets argument is neither a list like object (e.g. GRangesList),
           nor a set of files') 
  }
  
  
  # ----------------------------------------------------------------- #
  # checks whether the list argument contains only scoreMatrix objects
  if(all(unlist(lapply(targets, class)) == 'ScoreMatrix'))
    return(new("ScoreMatrixList",targets))
  
  # ----------------------------------------------------------------- #
  if(is.null(windows))
    stop("windows of class GRanges must be defined")
  
  # Given a list of RleList objects and a granges object, returns the scoreMatrix list Object
  if(list.ind & !all(unlist(lapply(targets, class)) %in% c('SimpleRleList', 'RleList','GRanges'))){
    stop('targets should be one of the following: an RleList, list of Rle, 
         GRangesList, a list of GRanges objects')
  }  
  
  if(!list.ind && all(file.exists(targets)) && is.null(type))
    stop('When providing a file, it is necessary to give the type of the file')
  
  # gets the names for the resulting list
  if(all(is.character(targets))){
    names = basename(targets)
  }else{
    names = names(targets)
  }
  
  sml = list()
  calc.ScoreMatrices <- function(i, targets=targets, windows=windows,
                                 strand.aware=strand.aware, weight.col=weight.col,
                                 is.noCovNA=is.noCovNA, type=type, rpm=rpm,
                                 unique=unique, extend=extend, param=param,
                                 bin.num,library.size=library.size){
    
    message(paste("working on", names(targets)[i]))
    
    if(is.null(bin.num) && all(width(windows) == unique(width(windows)))){
      
      ScoreMatrix(targets[[i]], windows=windows, 
                  strand.aware=strand.aware,
                  weight.col=weight.col, 
                  is.noCovNA=is.noCovNA,
                  type=type,
                  rpm=rpm,
                  unique=unique,
                  extend=extend,
                  param=param,
                  library.size=library.size[i])
    }else{
      if(is.null(bin.num))
        bin.num = 10
      
      ScoreMatrixBin(targets[[i]], windows=windows, 
                     bin.num=bin.num, 
                     bin.op=bin.op,
                     strand.aware=strand.aware,
                     weight.col=weight.col, 
                     is.noCovNA = is.noCovNA,
                     type=type,
                     rpm=rpm,
                     unique=unique,
                     extend=extend,
                     param=param,
                     library.size=library.size[i])
    }
  }    
  sml <- mclapply(1:length(targets), 
                  calc.ScoreMatrices,
                  targets=targets,
                  windows=windows,strand.aware=strand.aware,
                  weight.col=weight.col,is.noCovNA=is.noCovNA,type=type,
                  rpm=rpm,unique=unique,extend=extend,param=param,bin.num,
                  library.size=library.size,
                  mc.cores=cores)
  names(sml) = names
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
#' @return Shows the number of matrices and their sizes
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
#' Scales each ScoreMatrix in the ScoreMatrixList object, by rows and/or columns
#' 
#' @param sml a \code{ScoreMatrixList} object
#' @param columns a \code{columns} whether to scale the matrix by columns. 
#'               Set by default to FALSE
#' @param rows  a \code{rows} Whether to scale the matrix by rows. Set by default 
#'            to TRUE
#' @param scalefun a function object that takes as input a matrix and returns a 
#'        matrix.
#'  By default  the argument is set to the R scale function with center=TRUE and 
#'  scale=TRUE
#'
#' @usage scaleScoreMatrixList(sml, columns, rows, scalefun)
#' @return \code{ScoreMatrixList} object
#' @examples 
#' library(GenomicRanges)
#' data(cage)
#' data(cpgi)
#' data(promoters)
#'  
#' cage$tpm = NULL
#' targets = GRangesList(cage=cage, cpgi=cpgi)
#' sml = ScoreMatrixList(targets, promoters, bin.num=10, strand.aware=TRUE)
#' sml.scaled = scaleScoreMatrixList(sml, rows=TRUE)
#' sml.scaled
#' \donttest{
#' multiHeatMatrix(sml) 
#' }
#'
#' @docType methods
#' @rdname scaleScoreMatrixList
#' @export

setGeneric("scaleScoreMatrixList", 
           function(sml, 
                    columns=FALSE, rows=TRUE, 
                    scalefun=NULL) 
             standardGeneric("scaleScoreMatrixList") )

#' @aliases scaleScoreMatrixList,ScoreMatrixList-method
#' @rdname scaleScoreMatrixList
setMethod("scaleScoreMatrixList", signature("ScoreMatrixList"),
          function(sml, columns, rows, scalefun){
            
            
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
#' @examples
#' library(GenomicRanges)
#' target = GRanges(rep(c(1,2),each=7), 
#'                   IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5), 
#'                   weight = rep(c(1,2),each=7))
#'                  
#' windows1 = GRanges(rep(c(1,2),each=2), 
#'                     IRanges(rep(c(1,2), times=2), width=5), 
#'                     strand=c('-','+','-','+'))
#' windows2 = windows1[c(1,3)]
#' sml = as(list(ScoreMatrix(target, windows1),
#'                ScoreMatrix(target, windows2)), 'ScoreMatrixList')
#' sml
#' intersectScoreMatrixList(sml)
#'
#' @docType methods
#' @rdname intersectScoreMatrixList-methods
#' @export
setGeneric("intersectScoreMatrixList", 
           function(sml,reorder=FALSE)
             standardGeneric("intersectScoreMatrixList") )

#' @aliases intersectScoreMatrixList,ScoreMatrixList-method
#' @rdname intersectScoreMatrixList-methods
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
#' @return \code{ScoreMatrixList} object
#' 
#' @examples
#' library(GenomicRanges)
#' data(cage)
#' data(cpgi)
#' data(promoters)
#' 
#' cage$tpm = NULL
#' targets = GRangesList(cage=cage, cpgi=cpgi)
#' sml = ScoreMatrixList(targets, promoters, bin.num=10)
#' kmeans.clust = kmeans(sml$cage,3)
#'  
#' sml.ordered = orderBy(sml, kmeans.clust$cluster)
#' \donttest{
#' multiHeatMatrix(sml.ordered)
#' }
#' @docType methods
#' @rdname orderBy-methods
#' @export
setGeneric("orderBy", 
           function(sml,ord.vec)
             standardGeneric("orderBy") )

#' @aliases orderBy,ScoreMatrixList-method
#' @rdname orderBy-methods
setMethod("orderBy", signature("ScoreMatrixList"),
          function(sml, ord.vec){
            
            if(is.null(ord.vec))
              return(sml)
            
            if(!is.integer(ord.vec))
              stop('ord.vec needs to be of an integer vector')  
            
            sml = lapply(sml, function(x)x[ord.vec,])
            return (as(sml,'ScoreMatrixList'))
          }
)


# ---------------------------------------------------------------------------- #
#' @aliases binMatrix,ScoreMatrixList-method
#' @rdname binMatrix-methods
setMethod("binMatrix", signature("ScoreMatrixList"),
          function(x, bin.num=NULL, fun='mean'){
            
            if(is.null(bin.num))
              return(x)
            
            return(new("ScoreMatrixList", 
                       lapply(x, function(y)binMatrix(y, bin.num=bin.num, fun=fun))))
          }
)
