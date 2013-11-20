#######################################
# S3 functions
#######################################

# ---------------------------------------------------------------------------- #
### gets colors for a factor variable
getColors = function(n) {
  
  black = "#000000"
  c(black,hcl(h=seq(0,(n-2)/(n-1),length=n-1)*360,c=100,l=65,fixup=TRUE))
}

# ---------------------------------------------------------------------------- #
# removes ranges that fell of the rle object
# does not check for the correspondence of the chromosome names - always check before using this function
constrainRanges = function(target, windows){
	
  checkClass(target, c('SimpleRleList','RleList','CompressedRleList'))
	checkClass(windows, 'GRanges')
	
	mcols(windows)$X_rank = 1:length(windows)
	r.chr.len = elementLengths(target)
    constraint = GRanges(seqnames=names(r.chr.len),
                         IRanges(start=rep(1,length(r.chr.len)),
                                   end=as.numeric(r.chr.len)))
	# suppressWarnings is done becuause GenomicRanges function give warnings 
  #if you don't have the same seqnames in both objects
  win.list.chr = suppressWarnings(subsetByOverlaps(windows, 
                                                   constraint,
                                                   type = "within",
                                                   ignore.strand = TRUE))
	
	if(length(win.list.chr) == 0)
    stop('All windows fell have coordinates outside windows boundaries')
	return(win.list.chr)
}


# ---------------------------------------------------------------------------- #
# checkw whether the x object corresponds to the given class
checkClass = function(x, class.name, var.name = deparse(substitute(x))){

	fun.name = match.call(call=sys.call(sys.parent(n=1)))[[1]]
	if(!class(x) %in% class.name)
		stop(paste(fun.name,': ', 
               var.name, 
               ' is not of class: ', 
               paste(class.name, collapse=' '), 
               '\n', sep=''))
}


# ---------------------------------------------------------------------------- #
# given a big bam path reads the big wig file into a RleList
# to be used by ScoreMatrix:char,GRanges
readBam = function(target, windows, rpm=FALSE,
                   unique=FALSE, extend=0, param=NULL, ...){
 
  # generates the ScanBamParam object
  if(is.null(param)){
    param <- ScanBamParam(which=reduce(windows, ignore.strand=TRUE))
  }else{
    if(class(param) == 'ScanBamParam'){
      bamWhich(param) <- reduce(windows, ignore.strand=TRUE)
    }else{
      stop('param needs to be an object of clas ScanBamParam')
    }
  }
  
  # get the coverage vector for 
  # given locations
  alns <- granges(readGAlignmentsFromBam(target, param=param))# read alignments
  if(unique)
    alns = unique(alns)
  
  if(extend > 0)
    alns = resize(alns, width=extend, fix='start')
  if(extend < 0)
    stop('extend needs to be a positive integer')
  
  
  covs=coverage(alns)
  
  if(rpm){
    message('Normalizing to rpm ...')
    total = 1e6/sum(countBam(BamFile(target))$records)
    covs = covs * total
  }
  
  return(covs)
  
}

# ---------------------------------------------------------------------------- #
# given a big wig path reads the big wig file into a RleList
# to be used by ScoreMatrix:char,GRanges
readBigWig = function(target, windows=NULL, ...){
  
  
  if(class(windows) != 'GRanges')
      stop('windows argument needs to be a GRanges object')
  
  
  if(is.null(windows)){
    bw = import(target, asRangedData = FALSE)
  }else{
    bw = import(target, asRangedData = FALSE, which=windows)
  }
  if(length(bw) == 0)
    stop('There are no ranges selected')
  
  covs = coverage(bw, weight=bw$score)
  return(covs)
}



#######################################
# S4 functions
#######################################


# ---------------------------------------------------------------------------- #
#' Get base-pair score for bases in each window
#'
#' The funcion produces a base-pair resolution matrix of scores for given equal
#' width windows of interest. The returned matrix  can be used to 
#' draw meta profiles or heatmap of read coverage or wig track-like data.
#' The \code{windows} argument can be a predefined region around transcription start sites 
#' or other regions of interest that have equal lengths
#' The function removes all window that fall off the Rle object - 
#' have the start coordinate < 1 or end coordinate > length(Rle)
#' The function takes the intersection of names in the Rle and GRanges objects
#'
#' @param target \code{RleList} , \code{GRanges} or a BAM file
#'  to be overlapped with ranges in \code{windows}
#' @param windows \code{GRanges} object that contains the windows of interest. 
#'                It could be promoters, CpG islands, exons, introns. 
#'                However the sizes of windows have to be equal.
#' @param strand.aware If TRUE (default: FALSE), the strands of the
#'                   windows will be taken into account in the resulting
#'                    \code{ScoreMatrix}.
#'                     If the strand of a window is -, the values of the bins 
#'                     for that window will be reversed
#' @param weight.col if the object is \code{GRanges} object a numeric column
#'                 in meta data part can be used as weights. This is particularly
#'                useful when genomic regions have scores other than their
#'                coverage values, such as percent methylation, conservation
#'                scores, GC content, etc. 
#' @param is.noCovNA (Default:FALSE)
#'                  if TRUE,and if 'target' is a GRanges object with 'weight.col'
#'                   provided, the bases that are uncovered will be preserved as
#'                   NA in the returned object. This useful for situations where
#'                   you can not have coverage all over the genome, such as CpG methylation
#'                   values.
#' @param type if target is a character vector of file paths, then type designates
#'              the type of the corresponding files (bam or bigWig)
#' @param rpm boolean telling whether to normalize the coverage to per milion 
#'                    reads. FALSE by default.
#' @param unique boolean which tells the function to remove duplicated reads 
#'              based on chr, start, end and strand
#' @param extend numeric which tells the function to extend the reads to width=extend
#' @param param ScanBamParam object 
#' 
#' @return returns a \code{ScoreMatrix} object
#' @seealso \code{\link{ScoreMatrixBin}}
#' @examples
#' 
#' # When target is GRanges
#'  data(cage)
#'  data(promoters)
#'  scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
#'                                  weight.col="tpm")                                
#' # When target is RleList
#' covs = coverage(cage)
#' scores2 = ScoreMatrix(target=covs,windows=promoters,strand.aware=TRUE)    
#' 
#' # When target is a bam file
#'  bam.file = system.file('extdata/test.bam', package='genomation')
#'  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5))
#'  scores3 = ScoreMatrix(target=bam.file,windows=windows, type='bam') 
#'  
#' # when target is a bigWig file
#'  bw.file = system.file('tests/test.bw', package='rtracklayer')
#'  windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5))
#'  scores3 = ScoreMatrix(target=bam.file,windows=windows, type='bigWig') 
#'  
#' @docType methods
#' @rdname ScoreMatrix-methods           
#' @export
setGeneric("ScoreMatrix",
                    function(target,windows,
                             strand.aware=FALSE,
                             weight.col=NULL,
                             is.noCovNA=FALSE,
                             type='',
                             rpm=FALSE, 
                             unique=FALSE, 
                             extend=0,
                             param=NULL) 
                                standardGeneric("ScoreMatrix") )



#' @aliases ScoreMatrix,RleList,GRanges-method
#' @rdname ScoreMatrix-methods
setMethod("ScoreMatrix",signature("RleList","GRanges"),
          function(target,windows,strand.aware){
            
   #check if all windows are equal length
    if( length(unique(width(windows))) >1 ){
        stop("width of 'windows' are not equal, provide 'windows' with equal widths")
    }     
		
    # set a uniq id for the GRanges
    windows = constrainRanges(target, windows)
		
   
  	# fetches the windows and the scores
    chrs = intersect(names(target), as.character(unique(seqnames(windows))))
    myViews=Views(target[chrs],as(windows,"RangesList")[chrs]) # get subsets of RleList
    
    #  get a list of matrices from Views object
    #  operation below lists a matrix for each chromosome
    mat = lapply(myViews,function(x) t(viewApply(x,as.vector)) )
    #mat=as.matrix(myViews) # this might work as well - have to check which one is faster
    
    # combine the matrices from chromosomes 
    mat = do.call("rbind",mat)   
    
    # get the ranks of windows, when things are reorganized by as(...,"RangesList")
    r.list=split(mcols(windows)[,"X_rank"], as.factor(seqnames(windows))  )    
    ranks=do.call("c",r.list)
    
    # put original window order as rownames
    # this is very important to keep
    # if we have multiple RleLists for the same set of windows
    # we need to know which windows are removed from which set
    # so when we are doing something comparative (clustering windows
    # based on different histone marks) we only work with windows
    # that are covered by all histone marks
  	rownames(mat) = ranks  
    
    # if strand aware is TRUE, we need to flip the windows on the minus strand
  	if(strand.aware == TRUE){
      orig.rows=which(as.character(strand(windows)) == '-')
      mat[rownames(mat) %in% orig.rows,] = mat[rownames(mat) %in% 
                                               orig.rows, ncol(mat):1]
  	}

	# reorder matrix
  	mat = mat[order(ranks),] 
    
  return(new("ScoreMatrix",mat))
})



# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrix,GRanges,GRanges-method
#' @rdname ScoreMatrix-methods
setMethod("ScoreMatrix",signature("GRanges","GRanges"),
          function(target, windows, strand.aware, weight.col,is.noCovNA){
            
            #make coverage vector  from target
            if(is.null(weight.col)){
				      target.rle=coverage(target)
            }else{
              if(! weight.col %in% names(mcols(target)) ){
                  stop("provided column 'weight.col' does not exist in tartget\n")
              }
              if(is.noCovNA)
                { # adding 1 to figure out NA columns later
                  target.rle=coverage(target,weight=(mcols(target)[weight.col][,1]+1) )
                  mat=ScoreMatrix(target.rle,windows,strand.aware)
                  mat=mat-1 # substract 1
                  mat[mat<0]=NA # everything that are <0 are NA
                  return(mat)
                }
                target.rle=coverage(target,weight= weight.col ) 
                
            }
            

            # call ScoreMatrix function
            ScoreMatrix(target.rle,windows,strand.aware)
})

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrix,character,GRanges-method
#' @rdname ScoreMatrix-methods
setMethod("ScoreMatrix",signature("character","GRanges"),
          function(target,windows, strand.aware, type='', 
                   rpm=FALSE, unique=FALSE, extend=0, param=NULL){
            
            if(!file.exists(target)){
			      	stop("Indicated 'target' file does not exist\n")
            }
            
            fm = c('bam','bigWig')
            if(!type %in% fm)
              stop(paste(c('currently supported formats are', fm)))
            
            if(type == 'bam' & !grepl('bam$',target))
              warning('you have set type="bam", but the designated file does not have .bam extension')
            if(type == 'bigWig' & !grepl('bw$',target))
              warning('you have set type="bigWig", but the designated file does not have .bw extension')
            
            if(type == 'bam')
              covs = readBam(target, windows, rpm=rpm, unique=unique, 
                             extend=extend, param=param)
            if(type == 'bigWig')
              covs = readBigWig(target=target, windows=windows)            
            
            # get coverage vectors
            ScoreMatrix(covs,windows,strand.aware)
          })



# ---------------------------------------------------------------------------- #
#' Bins the columns of a matrix using a user provided function 
#'
#' @param x \code{ScoreMatrix} or a \code{ScoreMatrixList} object
#' @param bin.num \code{integer} number of bins in the final matrix
#' @param fun \code{character} vector or an anonymous function that will be used for binning
#'
#' @return \code{ScoreMatrix} or \code{ScoreMatrixList} object
#'
#' @examples
#' 
#' # binning the columns in a ScoreMatrix object
#' target = GRanges(rep(c(1,2),each=7), IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5), 
#' weight = rep(c(1,2),each=7), 
#' strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))
#' windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
#' strand=c('-','+','-','+'))
#' sm = ScoreMatrix(target, windows)
#' bin = binMatrix(sm, bin.num=2)
#'
#' @docType methods
#' @rdname binMatrix-methods
#' @export
setGeneric("binMatrix", 
              function(x, bin.num=NULL, fun='mean')
                standardGeneric("binMatrix") )

#' @aliases binMatrix,ScoreMatrix-method
#' @rdname binMatrix-methods
setMethod("binMatrix", signature("ScoreMatrix"),
			function(x, bin.num=NULL, fun='mean'){
		  
				if(is.null(bin.num))
					return(x)
					
				if(bin.num > ncol(x))
					stop("number of given bins is bigger 
                than the number of matrix columns")
		  
        if(is.character(fun))
				  fun = match.fun(fun)
        
				coord = binner(1, ncol(x), bin.num)
				bmat = mapply(function(a,b)apply(x[,a:b],1,fun), coord[1,], coord[2,])
										
				return(new("ScoreMatrix", bmat))
		 }
)

# ---------------------------------------------------------------------------- #
# show Methods
#' @rdname show-methods
#' 
setMethod("show", "ScoreMatrix",
          function(object){
            dims = dim(object)
            
            s=sprintf("  %s %d %d", "scoreMatrix with dims:", dims[1], dims[2])
            message(s)  
          })

# ---------------------------------------------------------------------------- #
#' Scales the values in the matrix by rows and/or columns
#'
#' @param mat \code{ScoreMatrix} object
#' @param columns \code{columns} whether to scale the matrix by columns. Set by default to FALSE.
#' @param rows  \code{rows} Whether to scale the matrix by rows. Set by default to TRUE
#' @param scalefun function object that takes as input a matrix and returns a matrix. 
#'         By default  the argument is set to (x - mean(x))/(max(x)-min(x)+1)
#'
#' @return \code{ScoreMatrix} object
#'
#' @examples
#' 
#' # scale the rows of a scoreMatrix object
#' target = GRanges(rep(c(1,2),each=7), IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5), 
#'           weight = rep(c(1,2),each=7), 
#'           strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))
#' windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
#'           strand=c('-','+','-','+'))
#' sm = ScoreMatrix(target, windows)
#' ssm = scaleScoreMatrix(sm, rows=TRUE)
#'
#' @docType methods
#' @rdname scaleScoreMatrix-methods
#' @export
setGeneric("scaleScoreMatrix", 
                function(mat, 
                         columns=FALSE, rows=TRUE, 
                         scalefun=NULL) 
                        standardGeneric("scaleScoreMatrix") )

#' @aliases scaleScoreMatrix,ScoreMatrix-method
#' @rdname scaleScoreMatrix-methods
setMethod("scaleScoreMatrix", signature("ScoreMatrix"),
          function(mat, columns, rows, scalefun){
            
            if(is.null(scalefun))
              scalefun = function(x)(x-mean(x))/(max(x)-min(x)+1)
            
            if(columns)
              mat = apply(mat, 2, scalefun)
            
            if(rows)
              mat = t(apply(mat,1,scalefun))
            
            return(new("ScoreMatrix", mat))
          }
)
