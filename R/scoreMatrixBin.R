.onUnload <- function (libpath) {
  library.dynam.unload("genomation", libpath)
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
#' @param target  \code{RleList}, \code{GRanges}, a BAM file or a bigWig file 
#'                object to be overlapped with ranges in \code{windows}
#' @param windows \code{GRanges} or \code{GRangesList} object that contains
#'                the windows of interest. It could be promoters, CpG islands, 
#'                exons, introns as GRanges object or GrangesList object representing
#'                exons of each transcript. Exons must be ordered by ascending rank
#'                by their position in transcript. The sizes of windows 
#'                does NOT have to be equal.
#' @param bin.num single \code{integer} value denoting how many bins there 
#'                should be for each window
#' @param bin.op bin operation that is either one of the following strings: 
#'                "max","min","mean","median","sum". The operation is applied on the 
#'                values in the bin. Defaults to "mean"
#' @param strand.aware If TRUE (default: FALSE), the strands of the windows will 
#'                     be taken into account in the resulting \code{scoreMatrix}. 
#'                     If the strand of a window is -, the values of the bins 
#'                     for that window will be reversed
#' @param weight.col if the object is \code{GRanges} object a numeric column
#'                in meta data part can be used as weights. This is particularly
#'                useful when genomic regions have scores other than their
#'                coverage values, such as percent methylation, conservation
#'                scores, GC content, etc. 
#' @param is.noCovNA (Default:FALSE)
#'                  if TRUE,and if 'target' is a GRanges object with 'weight.col'
#'                   provided, the bases that are uncovered will be preserved as
#'                   NA in the returned object. This useful for situations where
#'                   you can not have coverage all over the genome, such as CpG
#'                   methylation values.
#' @param type   (Default:"auto")
#'               if target is a character vector of file paths, then type designates 
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
#'                            is calculated using the Rsamtools idxstatsBam function:
#'                            sum(idxstatsBam(target)$mapped).
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
#' 
#' # Compute transcript coverage of a set of exons.
#' library(GenomicRanges)
#' bed.file = system.file("extdata/chr21.refseq.hg19.bed", 
#'                        package = "genomation")
#' gene.parts = readTranscriptFeatures(bed.file)
#' transcripts = split(gene.parts$exons, gene.parts$exons$name)
#' transcripts = transcripts[]
#' myMat3 = ScoreMatrixBin(target=cage, windows=transcripts[1:250], 
#'                     bin.num=10)
#' myMat3                     
#' @seealso \code{\link{ScoreMatrix}}
#' @docType methods
#' @useDynLib genomation
#' @importFrom Rcpp sourceCpp
#' @rdname ScoreMatrixBin-methods           
#' @export
setGeneric("ScoreMatrixBin",
           function(target,windows,
                    bin.num=10,bin.op="mean",
                    strand.aware=FALSE,
                    weight.col=NULL,is.noCovNA=FALSE,
                    type='auto',
                    rpm=FALSE,
                    unique=FALSE,
                    extend=0,
                    param=NULL,
                    bam.paired.end=FALSE,
                    library.size=NULL
           ) 
           standardGeneric("ScoreMatrixBin") )

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrixBin,RleList,GRanges-method
#' @rdname ScoreMatrixBin-methods
#' @usage \\S4method{ScoreMatrixBin}{RleList,GRanges}(target, windows, bin.num, bin.op,
#'                                                    strand.aware)
setMethod("ScoreMatrixBin",signature("RleList","GRanges"),
          function(target, windows, bin.num, bin.op, strand.aware){

            mcols(windows)$X_rank <- 1:length(windows)
            
            # check if windows have width > 1
            if( any(width(windows)==1) ){
              stop("provide 'windows' with widths greater than 1")
            } 
            
            # removes windows that fall of the chromosomes - window id is in values(windows)$X_rank 
            windows.len <- length(windows)
            windows <- constrainRanges(target, windows)
            if(length(windows)!=windows.len){
              warning(paste0(windows.len-length(windows),
                             " windows fall off the target"))
            }
            
            # checks whether some windows are shorter than the wanted window size
            wi <- IRanges::width(windows) < bin.num
            if(any(wi)){
              windows <- windows[!wi]
              if(length(windows) == 0)
                stop('all supplied windows have width < number of bins')
              warning('supplied GRanges object contains ranges of width < number of bins')
            }
            

            # fetches the windows and the scores
            chrs <- sort(intersect(names(target), as.character(unique(seqnames(windows)))))
            myViews <- Views(target[chrs],as(windows,"IntegerRangesList")[chrs]); # get subsets of coverage

 
            mat <- lapply(myViews,function(x) as.list((viewApply(x,as.vector,simplify=FALSE))))
            mat <- do.call("c",mat)
            
            functs = c("min",'mean','max','median','sum')
            if(!bin.op %in% functs)
              stop(paste(c('Supported binning functions are', functs,'\n')))
            if(bin.op =="mean"){
              mat_res <- listSliceMean(mat, bin.num)
            }else if(bin.op =="median"){
              mat_res <- listSliceMedian(mat, bin.num)
            }else if(bin.op =="sum"){
              mat_res <- listSliceSum(mat, bin.num)
            }else if(bin.op =="max"){
              mat_res <- listSliceMax(mat, bin.num)
            }else if(bin.op =="min"){
              mat_res <- listSliceMin(mat, bin.num)
            }
            
            # copied from scoreMatrix()
            # get the ranks of windows, when things are reorganized by as(...,"IntegerRangesList")
            r.list <- split(mcols(windows)[,"X_rank"], as.vector(seqnames(windows))  )
            r.list <- r.list[order(names(r.list))]
            ranks <- do.call("c",r.list)    
            rownames(mat_res) <- ranks
             
            # if strand aware is TRUE, we need to flip the windows on the minus strand
                  if(strand.aware == TRUE){
                    orig.rows=windows[strand(windows) == '-',]$X_rank
                    mat_res[rownames(mat_res) %in% orig.rows,] = mat_res[rownames(mat_res) %in% 
                                                                 orig.rows, ncol(mat_res):1]                                                                                                             
                  }
            # reorder matrix
            mat_res = mat_res[order(as.numeric(rownames(mat_res))),,drop=FALSE]

       new("ScoreMatrix", mat_res)
          })


# ---------------------------------------------------------------------------- #
#' @aliases  ScoreMatrixBin,GRanges,GRanges-method
#' @rdname ScoreMatrixBin-methods
#' @usage \\S4method{ScoreMatrixBin}{GRanges,GRanges}(target,windows,
#'                                                    bin.num,bin.op,
#'                                                    strand.aware,weight.col,
#'                                                    is.noCovNA)
setMethod("ScoreMatrixBin",signature("GRanges","GRanges"),
          function(target,windows,bin.num,bin.op,
                   strand.aware,weight.col,is.noCovNA){
            
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
#'                                                      bin.op='mean',strand.aware, 
#'                                                      is.noCovNA=FALSE, type='auto',
#'                                                      rpm, unique, extend, param,
#'                                                      bam.paired.end=FALSE, 
#'                                                      library.size=NULL)
setMethod("ScoreMatrixBin",signature("character","GRanges"),
          function(target, windows, bin.num=10, 
                   bin.op='mean', strand.aware, 
                   is.noCovNA=FALSE,
                   type='auto', rpm, unique, extend, param,
                   bam.paired.end=FALSE, library.size=NULL){
            
            if(!file.exists(target)){
              stop("Indicated 'target' file does not exist\n")
            }
            
            type = target.type(target, type)
            
            if( type=="bigWig" & rpm==TRUE)
              warning("rpm=TRUE is not supported for type='bigWig'")

            if(type == 'bam'){
              covs = readBam(target, windows, rpm=rpm, unique=unique, 
                             extend=extend, param=param,
                             paired.end=bam.paired.end, library.size=library.size)
              return(ScoreMatrixBin(covs,
                                    windows,
                                    bin.num=bin.num,
                                    bin.op=bin.op,
                                    strand.aware=strand.aware))
            }
            
            if(type == 'bigWig'){
              if(is.noCovNA==FALSE){
                covs = readBigWig(target=target, windows=windows)
                # get coverage vectors
                return(ScoreMatrixBin(covs,
                               windows,
                               bin.num=bin.num,
                               bin.op=bin.op,
                               strand.aware=strand.aware))
                
              }else{
                if(is.null(windows)){
                  target.gr = import(target)
                }else{
                  target.gr = import(target, which=windows)
                }
                if(length(target.gr) == 0)
                  stop('There are no ranges selected')
                
                # bigwig contains only one column for score
                weight.col=names(mcols(target.gr)) # name of a score vector
                
                return(ScoreMatrixBin(target.gr,
                               windows,
                               bin.num=bin.num,
                               bin.op=bin.op,
                               strand.aware=strand.aware,
                               weight.col=weight.col,
                               is.noCovNA=is.noCovNA))
              }
            }
          })

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrixBin,RleList,GRangesList-method
#' @rdname ScoreMatrixBin-methods
#' @usage  \\S4method{ScoreMatrixBin}{RleList,GRangesList}(target,windows,
#'                                                         bin.num, bin.op,
#'                                                         strand.aware)
setMethod("ScoreMatrixBin",signature("RleList","GRangesList"),
          function(target, windows, bin.num, bin.op, strand.aware){
            
            seqinfo(target) <- merge(seqinfo(target), seqinfo(windows))
            if (!isTRUEorFALSE(strand.aware))
              stop("'ignore.strand' must be TRUE or FALSE")
            
            if( length(names(windows)) != length(unique(names(windows))) )
              stop("Windows don't have unique names")
            
            #when the windows are unnamed, give the names
            if(is.null(names(windows)))
               names(windows) <- seq(1:length(windows))
            
            if(any(sum(width(windows))< bin.num))
              stop('Please remove transcripts that are shorter than bin.num')
            
            # Remove windows that fall of the chromosomes
            unlisted = unlist(windows)
            unlisted$'.trname' = rep(names(windows), elementNROWS(windows))
            if(is.unsorted(unlisted))
              unlisted = sort(unlisted)
            ex = genomation:::constrainRanges(target, 
                                              unlisted)
            # extracts the coverage vectors for individual exons and
            # concatenates them together
            ex_cvg = as(target[ex],'NumericList')
            dex.cvg = data.table(cvg = as.list(ex_cvg), id=ex$'.trname')
            dex.cvg = dex.cvg[,list(cvg=list(unlist(cvg, use.names=FALSE))),by=id]
            
            # reorders the coverage vectors to correspond to the original GRList
            dex.cvg = dex.cvg[match(unique(ex$'.trname'),dex.cvg$id),]
            
            if (strand.aware){
              dex.cvg$strand = unlist(runValue(strand(windows))[dex.cvg$id])
              dex.cvg[dex.cvg$strand == '-',cvg := list(list(rev(unlist(cvg)))), by=list(id)]
            }
            
            # constructs the bins for each transcript
            seqlen = dex.cvg[,length(unlist(cvg, use.names=TRUE)), by=id]$V1
            bins = IRangesList(lapply(seqlen,
                                       function(x)
                                         IRanges(breakInChunks(x, 
                                                               nchunk=bin.num))))
            names(bins) = dex.cvg$id
            rle = RleList(dex.cvg$cvg)
            names(rle) = dex.cvg$id
            
            # bins the signal
            my.vList = RleViewsList(
              lapply(names(rle),
                     function(seqname)
                       Views(rle[[seqname]], bins[[seqname]])))
            # Copied from genomation:::summarizeViewsRle
            # Calculate min/mean/max/median in each bin
            functs = c("min",'mean','max','median')
            if(!bin.op %in% functs)
              stop(paste(c('Supported binning functions are', functs,'\n')))
            if(bin.op=="min")
              sum.bins=unlist(viewMins(my.vList,na.rm=TRUE), use.names=FALSE)     
            if(bin.op=="max")
              sum.bins=unlist(viewMaxs(my.vList,na.rm=TRUE), use.names=FALSE)
            if(bin.op=="mean")
              sum.bins=unlist(viewMeans(my.vList,na.rm=TRUE), use.names=FALSE)
            if(bin.op=="median")
              sum.bins=unlist(viewApply(my.vList, 
                                        function(x) median(x, na.rm=TRUE), 
                                        simplify = TRUE), use.names=FALSE)
            
            mat = matrix( sum.bins, ncol=bin.num, byrow=TRUE )
            mat[is.nan(mat)] = NA
            new("ScoreMatrix",mat)
          })

# ---------------------------------------------------------------------------- #
#' @aliases  ScoreMatrixBin,GRanges,GRangesList-method
#' @rdname ScoreMatrixBin-methods
#' @usage \\S4method{ScoreMatrixBin}{GRanges,GRangesList}(target,windows,
#'                                                    bin.num,bin.op,
#'                                                    strand.aware,weight.col,
#'                                                    is.noCovNA)
setMethod("ScoreMatrixBin",signature("GRanges","GRangesList"),
          function(target,windows,bin.num,bin.op,
                   strand.aware,weight.col,is.noCovNA){
            
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
                target.rle=coverage(target,weight=weight.col)
              }
            }
            # call scoreMatrix function
            ScoreMatrixBin(target.rle,windows,bin.num,
                           bin.op,strand.aware)
          })

# ---------------------------------------------------------------------------- #
#' @aliases ScoreMatrixBin,character,GRangesList-method
#' @rdname ScoreMatrixBin-methods
#' @usage  \\S4method{ScoreMatrixBin}{character,GRangesList}(target, windows, bin.num=10,
#'                                                      bin.op='mean',strand.aware, 
#'                                                      weight.col=NULL,
#'                                                      is.noCovNA=FALSE, type='auto',
#'                                                      rpm, unique, extend, param,
#'                                                      bam.paired.end=FALSE, 
#'                                                      library.size=NULL)
setMethod("ScoreMatrixBin",signature("character","GRangesList"),
          function(target, windows, bin.num=10, 
                   bin.op='mean', strand.aware, 
                   weight.col=NULL, is.noCovNA=FALSE, type='auto', 
                   rpm, unique, extend, param,
                   bam.paired.end=FALSE, library.size=NULL){
            if(!file.exists(target)){
              stop("Indicated 'target' file does not exist\n")
            }
            type = target.type(target, type)
            if( type=="bigWig" & rpm==TRUE)
              warning("rpm=TRUE is not supported for type='bigWig'")
            
            if(type == 'bam'){
              covs = readBam(target, unlist(windows), 
                             rpm=rpm, unique=unique, 
                             extend=extend, param=param,
                             paired.end=bam.paired.end, library.size=library.size)
              return(ScoreMatrixBin(covs,
                                  windows,
                                  bin.num=bin.num,
                                  bin.op=bin.op,
                                  strand.aware=strand.aware))
            }
            if(type == 'bigWig'){
              if(is.noCovNA==FALSE){
                covs = readBigWig(target=target, windows=unlist(windows))
                # get coverage vectors
                return(ScoreMatrixBin(covs,
                                      windows,
                                      bin.num=bin.num,
                                      bin.op=bin.op,
                                      strand.aware=strand.aware))
                
              }else{
                if(is.null(windows)){
                  target.gr = import(target)
                }else{
                  target.gr = import(target, which=unlist(windows))
                }
                if(length(target.gr) == 0)
                  stop('There are no ranges selected')
                return(ScoreMatrixBin(target.gr,
                                      windows,
                                      bin.num=bin.num,
                                      bin.op=bin.op,
                                      strand.aware=strand.aware,
                                      weight.col=weight.col,
                                      is.noCovNA=is.noCovNA))
              }
            }
          })

