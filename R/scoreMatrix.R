#######################################
# S3 functions
#######################################

# ---------------------------------------------------------------------------- #
### gets colors for a factor variable
getColors = function(n) {
  
  black = "#000000"
  c(black,hcl(h=seq(0,(n-2)/(n-1),length=n-1)*360,c=100,l=65,fixup=TRUE))
}
### Extract file extension from file path
file.ext = function(x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}
### Check if a target file is in bam and bigWig formats
### by looking at the file extension
target.type = function(target, type=""){
  
  if(length(target)!=1){
    stop("target has to be a single character")
  }
  bw.exts = c("bw","bigWig","bigwig","BigWig", "BIGWIG", "BW")
  bam.exts = c("bam", "BAM", "Bam")
  if(type=="" | length(type)!=1 | nchar(type)==0 | is.null(type))
      stop(paste0('set argument type to "auto", "bam" or "bigWig"\n'))
  
  if(type=="auto"){
    # automatically recognize type of target files
    # by looking at the file extensions
    exts = file.ext(target)  
    if( any(!(exts %in% c(bam.exts, bw.exts))) )
      stop(paste0('currently supported formats are bam and bigWig\n'))
    if( exts %in% bw.exts)  return("bigWig") 
    if( exts %in% bam.exts) return("bam")
    
  }else if(type %in% c(bam.exts, bw.exts)){
    # if type is BAM or bigWig
    bam.grepl = paste0((paste0(bam.exts, collapse="$|")), "$")
    bw.grepl = paste0((paste0(bw.exts, collapse="$|")), "$")
    if(type %in% bam.exts){
      if(!grepl(bam.grepl, target)){
        warning(paste0('you have set type="',type,
                       '", but the designated file does not have ',
                       'an extension of BAM file (.bam)'))}
      return("bam")
    }
    if(type %in% bw.exts){
      if(!grepl(bw.grepl, target)){
        warning(paste0('you have set type="',type,
                       '", but the designated file does not have',
                       ' an extension of biWig file (.bw or .bigWig)')) }
      return("bigWig")
    }
  }
}

# ---------------------------------------------------------------------------- #
# removes ranges that fell of the rle object
# does not check for the correspondence of the chromosome names - always 
# check before using this function
constrainRanges = function(target, windows){
  
  checkClass(target, c('SimpleRleList','RleList','CompressedRleList'))
  checkClass(windows, 'GRanges')
  
  mcols(windows)$X_rank = 1:length(windows)
  r.chr.len = elementNROWS(target)
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
# check whether the x object corresponds to the given class
checkClass = function(x, class.name, var.name = deparse(substitute(x))){
  
  fun.name = match.call(call=sys.call(sys.parent(n=1)))[[1]]
  if(!class(x) %in% class.name)
    stop(paste(fun.name,': ', 
               var.name, 
               ' is not of class: ', 
               paste(class.name, collapse=' '), 
               '\n', sep=''))
}

galpTo2Ranges <- function(x)
{
  gr1 <- granges(first(x))
  gr2 <- granges(invertStrand(last(x)))
  ans <- split(c(gr1, gr2), rep.int(seq_along(x), 2L))
  names(ans) <- names(x)
  ans
}

# ---------------------------------------------------------------------------- #
# given a big bam path reads the big wig file into a RleList
# to be used by ScoreMatrix:char,GRanges
readBam = function(target, windows, rpm=FALSE,
                   unique=FALSE, extend=0, param=NULL, 
                   paired.end=FALSE, library.size=NULL, ...){

  # check the ScanBamParam object
  if(!is.null(param) & class(param)!='ScanBamParam')
    stop('param needs to be an object of class ScanBamParam')
  
  if(paired.end){
    
    flag <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE, 
                        isUnmappedQuery=FALSE, hasUnmappedMate=FALSE)
    # reads that map into window which stretches from start of the first window
    # to end of the last window
    if(is.null(param)){
      param <- ScanBamParam(which=reduce(windows, ignore.strand=TRUE), flag=flag)
    }else{
      bamWhich(param) <- reduce(windows, ignore.strand=TRUE)
    }
    alnp <- readGAlignmentPairs(target, param=param, use.names=TRUE)
    
    # remove rows that are duplicated when mates of pair map into two different windows
    alnp = alnp[!duplicated(names(alnp))]
    alns <- granges(alnp)
    
  }else{
    
    if(is.null(param)){
      param <- ScanBamParam(which=reduce(windows, ignore.strand=TRUE))
    }else{
      bamWhich(param) <- reduce(windows, ignore.strand=TRUE)
    }
    alns <- granges(readGAlignments(target, param=param, use.names=TRUE))
  }
  
  if(unique)
    alns = unique(alns)
  
  if(extend > 0)
    alns = resize(alns, width=extend, fix='start')
  if(extend < 0)
    stop('extend needs to be a positive integer')
  
  # get the coverage vector for given locations
  covs=coverage(alns)
  
  if(rpm){
    message('Normalizing to rpm ...')
    if(is.null(library.size)){
      total = 1e6/sum(idxStats(normalizePath(target))[3])
    }else{
      total = 1e6/library.size
    }
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
    bw = import(target)
  }else{
    bw = import(target, which=windows)
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
#' The function takes the intersection of names in the Rle and GRanges objects.
#' On Windows OS the function will give an error if the target is a file in .bigWig format.
#'
#' @param target \code{RleList} , \code{GRanges}, a BAM file or a BigWig
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
#'                   you can not have coverage all over the genome, such as CpG 
#'                   methylation values.
#' @param type  (Default:"auto")
#'              if target is a character vector of file paths, then type designates
#'              the type of the corresponding files (bam or bigWig).
#' @param rpm boolean telling whether to normalize the coverage to per milion 
#'                    reads. FALSE by default. See \code{library.size}.
#' @param unique boolean which tells the function to remove duplicated reads 
#'              based on chr, start, end and strand
#' @param extend numeric which tells the function to extend the reads to width=extend
#' @param param ScanBamParam object 
#' @param bam.paired.end boolean indicating whether given BAM file contains 
#'                       paired-end reads (default:FALSE).
#'                       Paired-reads will be treated as fragments.
#' @param library.size numeric indicating total number of mapped reads in a BAM file
#'                            (\code{rpm} has to be set to TRUE).
#'                            If is not given (default: NULL) then library size 
#'                            is calculated using a Samtools idxstats like function:
#'                            sum(idxStats(target)$mapped).
#' 
#' @note
#' We assume that a paired-end BAM file contains reads with unique ids and we remove 
#' both mates of reads if they are repeated. Due to the fact that \code{ScoreMatrix} 
#' uses the GenomicAlignments:readGAlignmentPairs function to read paired-end BAM files
#' a duplication of reads occurs when mates of one pair map into two different windows.
#' 
#' Strands of reads in a paired-end BAM are inferred depending on strand of 
#' first alignment from the pair. This is a default setting in the 
#' GenomicAlignments:readGAlignmentPairs function (see a strandMode argument). 
#' This mode should be used when the paired-end data was generated using 
#' one of the following stranded protocols: 
#' Directional Illumina (Ligation), Standard SOLiD.
#' 
#' @return returns a \code{ScoreMatrix} object
#' @seealso \code{\link{ScoreMatrixBin}}
#' 
#' @examples
#
#' # When target is GRanges
#' data(cage)
#' data(promoters)
#' scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
#'                                 weight.col="tpm")                                
#'                                  
#' 
#' # When target is RleList
#' library(GenomicRanges)
#' covs = coverage(cage)
#' scores2 = ScoreMatrix(target=covs,windows=promoters,strand.aware=TRUE)   
#' scores2 
#'
#' # When target is a bam file
#' bam.file = system.file('unitTests/test.bam', package='genomation')
#' windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5))
#' scores3 = ScoreMatrix(target=bam.file,windows=windows, type='bam') 
#' scores3
#' 
#' @useDynLib genomation
#' @importFrom Rcpp sourceCpp
#' @docType methods
#' @rdname ScoreMatrix-methods           
#' @export
setGeneric("ScoreMatrix",
           function(target,windows,
                    strand.aware=FALSE,
                    weight.col=NULL,
                    is.noCovNA=FALSE,
                    type="auto",
                    rpm=FALSE, 
                    unique=FALSE, 
                    extend=0,
                    param=NULL,
                    bam.paired.end=FALSE,
                    library.size=NULL) 
             standardGeneric("ScoreMatrix") )


#' @aliases ScoreMatrix,RleList,GRanges-method
#' @rdname ScoreMatrix-methods
#' @usage  \\S4method{ScoreMatrix}{RleList,GRanges}(target,windows,strand.aware)
setMethod("ScoreMatrix",signature("RleList","GRanges"),
          function(target,windows,strand.aware){
            
            #check if all windows are equal length
            if( length(unique(width(windows))) >1 ){
              stop("width of 'windows' are not equal, provide 'windows' with equal widths")
            }     
            
            # set a uniq id for the GRanges
            windows.len=length(windows)
            windows = constrainRanges(target, windows)
            if(length(windows)!=windows.len){
              warning(paste0(windows.len-length(windows),
                             " windows fall off the target"))
            }
            
            # fetches the windows and the scores
            chrs = sort(intersect(names(target), as.character(unique(seqnames(windows)))))
            myViews=Views(target[chrs],as(windows,"RangesList")[chrs]) # get subsets of RleList
            
            #  get a list of matrices from Views object
            #  operation below lists a matrix for each chromosome
            mat = lapply(myViews,function(x) t(viewApply(x,as.vector)) )
            
            # combine the matrices from chromosomes 
            mat = do.call("rbind",mat)   
            
            # get the ranks of windows, when things are reorganized by as(...,"RangesList")
            r.list=split(mcols(windows)[,"X_rank"], as.vector(seqnames(windows))  )
            r.list=r.list[order(names(r.list))]
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
              orig.rows=windows[strand(windows) == '-',]$X_rank
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
#' @usage \\S4method{ScoreMatrix}{GRanges,GRanges}(target, windows, strand.aware,
#'                                                 weight.col, is.noCovNA)
setMethod("ScoreMatrix",signature("GRanges","GRanges"),
          function(target, windows, strand.aware, weight.col, is.noCovNA){
            
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
#' @usage \\S4method{ScoreMatrix}{character,GRanges}(target, windows, strand.aware, 
#'                                                   type="auto", rpm=FALSE,
#'                                                   unique=FALSE, extend=0, param=NULL, 
#'                                                   bam.paired.end=FALSE,
#'                                                   library.size=NULL)
setMethod("ScoreMatrix",signature("character","GRanges"),
          function(target,windows, strand.aware, type='auto', 
                   rpm=FALSE, unique=FALSE, extend=0, 
                   param=NULL, bam.paired.end=FALSE,
                   library.size=NULL){
            
            if(!file.exists(target)){
              stop("Indicated 'target' file does not exist\n")
            }
            
            type = target.type(target, type)
            
            if( type=="bigWig" & rpm==TRUE)
              warning("rpm=TRUE is not supported for type='bigWig'")
            
            if(type == 'bam')
              covs = readBam(target, windows, rpm=rpm, unique=unique, 
                             extend=extend, param=param, 
                             paired.end=bam.paired.end, library.size)
            if(type == 'bigWig')
              covs = readBigWig(target=target, windows=windows)            
            
            #get coverage vectors
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
#' library(GenomicRanges)
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
#' @return Shows the dimension of the ScoreMatrix
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
#' library(GenomicRanges)
#' target = GRanges(rep(c(1,2),each=7), IRanges(rep(c(1,1,2,3,7,8,9), times=2), width=5), 
#'          weight = rep(c(1,2),each=7), 
#'          strand=c('-', '-', '-', '-', '+', '-', '+', '-', '-', '-', '-', '-', '-', '+'))
#' windows = GRanges(rep(c(1,2),each=2), IRanges(rep(c(1,2), times=2), width=5), 
#'            strand=c('-','+','-','+'))
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
