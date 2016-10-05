#######################################
# S3 functions
#######################################
# ---------------------------------------------------------------------------- #
# given a vector and length smooths the vector to a given size
# the function is not safe - check for the window length before
binner=function(start,end,nbins){
  
  if(! is.numeric(start))
    stop('start needs to be class numeric')
  if(! is.numeric(end))
    stop('end needs to be class numeric')
  if(! is.numeric(nbins))
    stop('nbins needs to be class numeric')
  
  x = unique(seq(from = start, to = end,length.out=nbins + 1 ) )
  my.start = ceiling(x)[-length(x)]
  my.end = floor(x)[-1]
  
  return( t(cbind(my.start, my.end) )  )
}

# ---------------------------------------------------------------------------- #
# given a target Rle and windows gets the views to be used for binning
getViewsBin = function(target, windows, bin.num){
  
  #get coordinates of bins in each window
  coord = matrix(mapply(binner, 
                        IRanges::start(windows),
                        IRanges::end(windows), 
                        bin.num, SIMPLIFY=TRUE), 
                 ncol=2, byrow=TRUE)
  
  # make GRanges object for the bins
  # subtract 1 so next start pos is not identical to  current end pos
  # keep window rank from original "windows" GRanges object
  subWins = GRanges(seqnames=rep(as.character(seqnames(windows)),each=bin.num),
                    IRanges(start=coord[,1],end=coord[,2])) 
  IRanges::values(subWins)$X_rank = rep(IRanges::values(windows)$X_rank, each=bin.num)
  
  # convert sub-windows to RangesList to be fed into coverage()
  win.list=as(subWins, "RangesList")
  win.list = win.list[sapply(win.list, length) > 0] # remove chr with no views on
  
  #check if there are common chromsomes
  chrs  = intersect(names(win.list), names(target))
  if(length(chrs)==0)
    stop("There are no common chromosomes/spaces to do overlap")
  
  # get views on your windows
  my.vList = Views(target[chrs], win.list[chrs] )
  my.vList = lapply(chrs, 
                    function(x){
                      v = my.vList[[x]]
                      names(v) = IRanges::values(subWins)$X_rank[as.character(seqnames(subWins)) == x]
                      return(v)})
  names(my.vList) = chrs
  return(my.vList)
}


# ---------------------------------------------------------------------------- #
# applies the summary function for the views to bin the objects - for standard Rle and returns a matrix object
summarizeViewsRle = function(my.vList, windows, bin.op, bin.num, strand.aware){
  
  # chop windows to bins
  functs = c("min",'mean','max','median')
  if(!bin.op %in% functs)
    stop(paste(c('Supported binning functions are', functs,'\n')))
  if(bin.op=="min")
    sum.bins=unlist(IRanges::lapply(my.vList, function(x) 
      IRanges::viewMins(x,na.rm=TRUE) ),use.names=FALSE)      
  
  if(bin.op=="max")
    sum.bins=unlist(IRanges::lapply(my.vList, function(x) 
      IRanges::viewMaxs(x,na.rm=TRUE) ),use.names=FALSE)      
  if(bin.op=="mean")
    sum.bins=unlist(IRanges::lapply(my.vList, function(x) 
      IRanges::viewMeans(x,na.rm=TRUE) ),use.names=FALSE)    
  
  if(bin.op=="median")
    sum.bins=unlist(IRanges::lapply(my.vList, function(x) 
      viewApply(x, function(x) median(as.numeric(x),na.rm=TRUE)  )), use.names=FALSE) 
  
  mat=matrix( sum.bins, ncol=bin.num,byrow=TRUE)
  mat[is.nan(mat)]=NA
  rownames(mat) = unlist(IRanges::lapply(my.vList, names), use.names=FALSE)[seq(1, length(mat), bin.num)]
  if(strand.aware){
    orig.rows=windows[strand(windows) == '-',]$X_rank
    mat[rownames(mat) %in% orig.rows,] = mat[rownames(mat) %in% orig.rows, ncol(mat):1]
  }
  mat = mat[order(as.numeric(rownames(mat))),]
  return(mat)
  
}


#######################################
# S4 functions
#######################################

# ---------------------------------------------------------------------------- #
#' Get bin score for bins on each window
#'
#' The function first bins each window to equal number of bins, and calculates
#' the a summary matrix for scores of each bin (currently, mean, max and min supported)
#' A scoreMatrix object can be used to draw average profiles or heatmap of read coverage
#' or wig track-like data. \code{windows} can be a predefined region such as CpG islands,
#' gene bodies, transcripts or CDS (coding sequences) that are not necessarily equi-width. 
#' Each window will be chopped to equal number of bins based on \code{bin.num} option.
#'
#' @param target  \code{RleList},  \code{GRanges}, BAM file or a bigWig file 
#'               object to be overlapped with ranges in \code{windows}
#' @param windows \code{GRanges} object that contains the windows of interest. 
#'                It could be promoters, CpG islands, exons, introns. However, 
#'                the sizes of windows does NOT have to be equal.
#' @param bin.num single \code{integer} value denoting how many bins there 
#'                should be for each window
#' @param bin.op bin operation that is either one of the following strings: 
#'              "max","min","mean". The operation is applied on the 
#'              values in the bin. Defaults to "mean"
#' @param strand.aware If TRUE (default: FALSE), the strands of the windows will 
#'                     be taken into account in the resulting \code{scoreMatrix}. 
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
#'                    methylation values.
#' @param type if target is a character vector of file paths, then type designates 
#'               the type of the corresponding files (bam or bigWig)
#' @param rpm boolean telling whether to normalize the coverage to per milion reads. 
#'            FALSE by default. See \code{library.size}.
#' @param unique boolean which tells the function to remove duplicated reads 
#'               based on chr, start, end and strand
#' @param extend numeric which tells the function to extend the reads to width=extend
#' @param param ScanBamParam object 
#' @param bam.paired.end boolean indicating whether given BAM file contains paired-end 
#'                       reads (default:FALSE). Paired-reads will be treated as 
#'                       fragments.             
#' @param library.size numeric indicating total number of mapped reads in a BAM file
#'                            (\code{rpm} has to be set to TRUE).
#'                            If is not given (default: NULL) then library size 
#'                            is calculated using a Samtools idxstats like function:
#'                            sum(idxStats(target)$mapped).
#'                                          
#' @return returns a \code{scoreMatrix} object
#' 
#' @examples
#' data(cage)
#' data(cpgi)
#' data(promoters)
#' myMat=ScoreMatrixBin(target=cage,
#'                       windows=cpgi,bin.num=10,bin.op="mean",weight.col="tpm")
#' \donttest{
#' plot(colMeans(myMat,na.rm=TRUE),type="l")
#' }
#'   
#' myMat2=ScoreMatrixBin(target=cage,
#'                        windows=promoters,bin.num=10,bin.op="mean",
#'                        weight.col="tpm",strand.aware=TRUE)
#' \donttest{
#' plot(colMeans(myMat2,na.rm=TRUE),type="l")
#' }
#' @seealso \code{\link{ScoreMatrix}}
#' @docType methods
#' @rdname ScoreMatrixBin-methods           
#' @export
setGeneric("ScoreMatrixBin",
           function(target,windows,
                    bin.num=10,bin.op="mean",
                    strand.aware=FALSE,
                    weight.col=NULL,is.noCovNA=FALSE,
                    type='',
                    rpm=FALSE,
                    unique=FALSE,
                    extend=0,
                    param=NULL,
                    bam.paired.end=FALSE,
                    library.size=NULL
           ) 
           standardGeneric("ScoreMatrixBin") )


# ---------------------------------------------------------------------------- #

# TODO: check if exons within the same granges in grl are form the same chromosome
# two exons from the same gene shouldne be on two different chromosomes

#' @aliases ScoreMatrixBin,RleList,GRangesList-method
#' @rdname ScoreMatrixBin-methods
#' @usage  \\S4method{ScoreMatrixBin}{RleList,GRangesList}(target,windows,strand.aware)
setMethod("ScoreMatrixBin",signature("RleList","GRangesList"),
          function(target, windows, bin.num, bin.op, strand.aware){
            
            # RleList target, GRangesList windows (exons)
            
            seqinfo(target) <- merge(seqinfo(target), seqinfo(windows))
            if (!isTRUEorFALSE(strand.aware))
              stop(msg("'ignore.strand' must be TRUE or FALSE"))  
            
            # Remove windows that fall of the chromosomes
            # Window ids are in values(windows)$X_rank 
            wind.nams = names(windows)
            wind.eltNROWS = elementNROWS(windows)
            
            unlisted <- unlist(windows)
            extra.col = rep(wind.nams, wind.eltNROWS)
            unlisted$transcript.name <- extra.col
            
            ex = genomation:::constrainRanges(target, 
                                              unlisted )
            
            tbl = table(ex$transcript.name) # TODO: is it fast enough?
            n.elements.ingrl = as.numeric(tbl) 
            names.elements.ingrl = names(tbl)
            windows = relist(
              ex, 
              # PartitioningByWidth indicates only structure of
              # new GRangesList after removing windows out of chrs
              # content (start, end, strand) doesnt matter for 
              # the relist function
              PartitioningByWidth( n.elements.ingrl )
            )
            names(windows) <- names.elements.ingrl
            
            # Checks whether some windows are shorter than the wanted window size
            wi = IRanges::width(ex) < bin.num
            if(any(wi)){
              windows = windows[!wi]
              if(length(windows) == 0)
                stop('all supplied windows have width < number of bins')
              warning('supplied GRanges object contains ranges of width < number of bins')
            }
            
            # Copied from GenomicFeatures:coverageByTranscript
            ## We could simply do 'uex <- unique(ex)' here but we're going to need
            ## 'sm' and 'is_unique' later to compute the "reverse index" so we compute
            ## them now and use them to extract the unique exons. That way we hash
            ## 'ex' only once (the expensive operation).
            sm <- selfmatch(ex)  # uses a hash table internally
            is_unique <- sm == seq_along(sm)
            uex2ex <- which(is_unique)  # index of unique exons
            uex <- ex[uex2ex]  # unique exons
            
            ## 2) Compute coverage for each unique exon ('uex_cvg').
            # coverage for every window regardless from which transcript they come from
            cvg <- target
            uex_cvg <- cvg[uex]  # parallel to 'uex'
            
            # This step is not needed because
            # they are already reverted by genomation in summarizeViewsRle function
            ## 3) Flip coverage for exons on minus strand.
            ## It feels like this is not as fast as it could be (the bottleneck being
            ## subsetting an Rle object which needs to be revisited at some point).
            #uex_cvg <- revElements(uex_cvg, strand(uex) == "-") 

            ## 4) Compute coverage by original exon ('ex_cvg').
            ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]  # reverse index
            stopifnot(identical(ex2uex[uex2ex], seq_along(uex2ex)))  # sanity
            stopifnot(identical(ex2uex[sm], ex2uex))  # sanity
            stopifnot(all(uex[ex2uex] == ex))  # sanity
            
            ex_cvg <- uex_cvg[ex2uex]  # parallel go 'ex'
            
            ## 5) Compute coverage of each transcript by concatenating coverage of its
            ##    exons.
            ans <- IRanges:::regroupBySupergroup(ex_cvg, windows)
            
            ## 6) Propagate 'mcols(transcripts)'.
            mcols(ans) <- mcols(windows) 
            #names(ans) <- names(windows) # if not mcols then names works.
            
            # gets the view list
            rlist = lapply(1:length(windows),
                           function(i){
                             my.vList.i = genomation:::getViewsBin(target, windows[[i]], bin.num)
                             #  class(my.vList.i) list
                             genomation:::summarizeViewsRle(my.vList.i, windows[[i]],
                                                            bin.op, bin.num, strand.aware)
                           })
            mat = do.call(rbind, rlist)
            new("ScoreMatrix",mat)
          })

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrixBin,character,GRangesList-method
#' @rdname ScoreMatrixBin-methods
#' @usage  \\S4method{ScoreMatrixBin}{character,GRangesList}(target,windows,strand.aware)
setMethod("ScoreMatrixBin",signature("character","GRangesList"),
          function(target, windows, bin.num=10, 
                   bin.op='mean', strand.aware, 
                   type, rpm, unique, extend, param,
                   bam.paired.end=FALSE, library.size=NULL){
            
            print("COSSS")
            
            if(!file.exists(target)){
              stop("Indicated 'target' file does not exist\n")
            }
            
            fmbw = c('bigWig','bw','bigwig','BigWig')
            if(!type %in% c('bam', fmbw)){
              if(type==""){
                stop(paste0('set argument type to "bam" or "bigWig"\n'))
              }
              stop('currently supported formats are bam and bigWig\n')
            }
            
            if(type == 'bam' & !grepl('bam$',target))
              warning('you have set type="bam", but the designated file does not have .bam extension')
            if(type == 'bigWig' & !grepl('bw$|bigWig$|bigwig$|BigWig$',target))
              warning('you have set type="bigWig", but the designated file does not have .bw extension')
            
            print("AAAAAAAAAAAAAAAAA")
            if(type == 'bam')
              covs = readBam(target, windows, rpm=rpm, unique=unique, 
                             extend=extend, param=param,
                             paired.end=bam.paired.end, library.size=library.size)
            if(type %in% fmbw)
              covs = readBigWig(target=target, windows=windows)        
            
            print("2. COSSSSSS")
            # get coverage vectors
            ScoreMatrixBin(covs,
                           windows,
                           bin.num=bin.num,
                           bin.op=bin.op,
                           strand.aware=strand.aware)
          })

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrixBin,RleList,GRanges-method
#' @rdname ScoreMatrixBin-methods
#' @usage \\S4method{ScoreMatrixBin}{RleList,GRanges}(target, windows, bin.num, bin.op,
#'                                                    strand.aware)
setMethod("ScoreMatrixBin",signature("RleList","GRanges"),
          function(target, windows, bin.num, bin.op, strand.aware){
            
            # removes windows that fall of the chromosomes - window id is in values(windows)$X_rank 
            windows = constrainRanges(target, windows)
            
            # checks whether some windows are shorter than the wanted window size
            wi = IRanges::width(windows) < bin.num
            if(any(wi)){
              windows = windows[!wi]
              if(length(windows) == 0)
                stop('all supplied windows have width < number of bins')
              warning('supplied GRanges object contains ranges of width < number of bins')
            }
            
            # gets the view list
            my.vList = getViewsBin(target, windows, bin.num)
            
            # summarize views with the given function
            mat = summarizeViewsRle(my.vList, windows, bin.op, bin.num, strand.aware)
            new("ScoreMatrix",mat)
          })


# ---------------------------------------------------------------------------- #
#' @aliases  ScoreMatrixBin,GRanges,GRanges-method
#' @rdname ScoreMatrixBin-methods
#' @usage \\S4method{ScoreMatrixBin}{GRanges,GRanges}(target,windows,bin.num,bin.op,
#'                                                    strand.aware,weight.col,is.noCovNA)
setMethod("ScoreMatrixBin",signature("GRanges","GRanges"),
          function(target,windows,bin.num,bin.op,strand.aware,weight.col,is.noCovNA){
            
            #make coverage vector  from target
            if(is.null(weight.col)){
              target.rle=coverage(target)
            }else{
              if(! weight.col %in% names(mcols(target)) ){
                stop("provided column 'weight.col' does not exist in target\n")
              }
              if(is.noCovNA)
              {  
                # adding 1 to figure out NA columns later
                target.rle=coverage(target,weight=(mcols(target)[weight.col][,1]+1) )
                
                # figure out which ones are real 0 score
                # which ones has no coverage
                # put NAs for no coverage bases
                runValue(target.rle)=runValue(target.rle)-1
                runValue(target.rle)[runValue(target.rle)<0]=NA
                
              }else{
                # get coverage with weights
                target.rle=coverage(target,weight= weight.col ) 
              }
            }
            
            # call scoreMatrix function
            ScoreMatrixBin(target.rle,windows,bin.num,
                           bin.op,strand.aware)
            
          })

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrixBin,character,GRanges-method
#' @rdname ScoreMatrixBin-methods
#' @usage \\S4method{ScoreMatrixBin}{character,GRanges}(target, windows, bin.num=10,
#'                                                      bin.op='mean',strand.aware, type,
#'                                                      rpm, unique, extend, param,
#'                                                      bam.paired.end=FALSE, 
#'                                                      library.size=NULL)
setMethod("ScoreMatrixBin",signature("character","GRanges"),
          function(target, windows, bin.num=10, 
                   bin.op='mean', strand.aware, 
                   type, rpm, unique, extend, param,
                   bam.paired.end=FALSE, library.size=NULL){
            
            if(!file.exists(target)){
              stop("Indicated 'target' file does not exist\n")
            }
            
            fmbw = c('bigWig','bw','bigwig','BigWig')
            if(!type %in% c('bam', fmbw)){
              if(type==""){
                stop(paste0('set argument type to "bam" or "bigWig"\n'))
              }
              stop('currently supported formats are bam and bigWig\n')
            }
            
            if(type == 'bam' & !grepl('bam$',target))
              warning('you have set type="bam", but the designated file does not have .bam extension')
            if(type == 'bigWig' & !grepl('bw$|bigWig$|bigwig$|BigWig$',target))
              warning('you have set type="bigWig", but the designated file does not have .bw extension')
            
            if(type == 'bam')
              covs = readBam(target, windows, rpm=rpm, unique=unique, 
                             extend=extend, param=param,
                             paired.end=bam.paired.end, library.size=library.size)
            if(type %in% fmbw)
              covs = readBigWig(target=target, windows=windows)        
            
            # get coverage vectors
            ScoreMatrixBin(covs,
                           windows,
                           bin.num=bin.num,
                           bin.op=bin.op,
                           strand.aware=strand.aware)
          })
