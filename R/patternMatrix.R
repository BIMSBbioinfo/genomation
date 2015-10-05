
# Notes:
# windows can have "N" but then score is set to 0 by
# patterns cannot contain "N"

# quantile(p@.Data, probs = c(0, 0.25, 0.5, 0.75, 0.9))

# BSgenome:matchPattern

# Questions:
# if windows contain only N then return matrix with 0s or nothing (right now nothing)

#For PWMs, we can either return occurence if the user provides a score cutoff, 
#or return scores themselves, I think the other function, motifScanHits() is capable of doing both
# TODO:  some basic benchmarks for 1000,5000,10000 windows



#######################################
# S3 functions
#######################################

# adapted from seqPattern::motifScanHits
# it returns granges instead of start(granges)
.scan.sequence.with.pwm <- function(pwm, seq, minScore){
  
  nc <- nchar(minScore)
  if (substr(minScore, nc, nc) == "%"){
    perc.threshold <- substr(minScore, 1L, nc-1L)
    min.score <- minScore(pwm)
    max.score <- maxScore(pwm)
    score.threshold = min.score + as.double(perc.threshold)/100 *
      (max.score-min.score)
  }else{
    score.threshold <- minScore
  }
  
  # make all PWM scores positive to avoid matching the Ns which are assigned score 0
  pwm.corr <- pwm + abs(min(pwm))
  score.threshold <- score.threshold + ncol(pwm) * abs(min(pwm))
  
  pwm.match <- matchPWM(pwm = pwm.corr, subject = seq, min.score = score.threshold)
  return(pwm.match)
}



.patternDNAStringSet_windowsDNAStringSet2matrix = function(pattern, windows, cores){


	   # Improve it for ScoreMatrixList
	   if(cores>1){
              # seqPattern::getPatternOccurrenceList uses parallel::mclapply for parallelization
              pat = getPatternOccurrenceList(windows, pattern,
                                             useMulticore=TRUE, nrCores=cores)
            }else{
              pat = getPatternOccurrenceList(windows, pattern)
            }
            
              pat = pat[[1]] #first 'pattern' -> for ScoreMatrix, if more then sml 
              
              if(is.null(pat)){
		              stop("Windows contain letters not in [ACGT], e.g. contain only 'N'")
              }
              
              gr <- GRanges(seqnames = "chrN",
                            ranges = IRanges(start=pat$position, width=length(pattern[[1]])),
                            strand = "*",
                            sequence = pat$sequence,
                            seqinfo = Seqinfo(seqnames="chrN",
                                              seqlengths=length(windows[[1]])))
              grl <- split(gr, gr$sequence)

              # calculate coverage of occurence of each pattern, 1-present, 0-absent
              grl.reduce <- reduce(grl)
              ma <- lapply(1:length(windows), function(x) as.vector(coverage(grl.reduce[x])[[1]]) )
              # convert list to matrix
              mat <- matrix(unlist(ma), ncol = length(ma[[1]]), byrow = TRUE)
              return(mat)

}

.patternPWM_windowsDNAStringSet2matrix = function(pattern, windows, 
                                                  min.score, asPercentage){

	    if(is.null(min.score)){
              mat <- motifScanScores(regionsSeq = windows,
                                             motifPWM = pattern,
                                             asPercentage = asPercentage)
            }else{
                            
              pwm.match <- lapply(1:length(windows), function(x){
                .scan.sequence.with.pwm(pwm = pattern, seq = toString(windows[x]), minScore=min.score) 
              })

              #ma <- lapply(1:length(windows), function(x) as.vector(coverage(pwm.match[[x]])) )
              # calculate coverage of occurence for each pattern, 1-present, 0-absent
              ma <- lapply(1:length(windows), function(x) as.vector(coverage(reduce(pwm.match[[x]]))) )
              
              # convert list to matrix
              mat <- matrix(unlist(ma), ncol = length(ma[[1]]), byrow = TRUE)
              
            }
            
            return(mat)
}



#######################################
# S4 functions
#######################################

# ---------------------------------------------------------------------------- #
#' Get scores that correspond to k-mer occurence for bases in each window
#'
#' we are looking for occurrence, we aim to get a scoreMatrix with 1s and 0s for k-mers.
#' The function first look for and calculates
#' The function firsts bins each window to equal number of bins, and calculates
#' the a summary matrix for scores of each bin (currently, mean, max and min supported)
#' A patternMatrix object can be used to ..
#' \code{windows} can be a predefined region such as CpG islands or gene bodies that have the same width.
#'
#' @param pattern matrix, a PWM matrix. It has one row for each nucleotide ("A","C","G" and "T")
#' @param windows \code{GRanges} object 
#' @param genome \code{BSgenome} object
#' @param min.score numeric or string indicating minimum score for counting a match. 
#'                  Can be given as a character string containing
#'                  a percentage of the highest possible score or as a single number 
#'                  (as default "80%"). If min.score is not provided or is NULL
#'                  then \code{patternMatrix} returns scores.
#' @param asPercentage boolean telling whether scores represent percentage of the maximal 
#'                     motif PWM score (default: TRUE) or raw scores (FALSE).
#' @param winsorize                  
#' @param cores the number of cores to use (default: 1)   
#'
#' @details
#' we are looking for occurrence, we aim to get a scoreMatrix with 1s and 0s for k-mers.
#' If one pattern is provided then created is a ScoreMatrix object, otherwise a ScoreMatrixList 
#' object.             
#' if \code{min.score} argument is NULL then returns scores (seqPattern:motifScanScores)
#' otherwise 0 or 1 (function adapted from seqPattern::motifScanHits)
#' 
#' @return returns a \code{scoreMatrix} object
#' 
#' @seealso \code{\link{ScoreMatrix}}
#' @docType methods
#' @rdname patternMatrix-methods           
#' @export
#' 
setGeneric(
  name="patternMatrix",
  def=function(pattern, windows, 
               genome=NULL, min.score = NULL, 
               cores=1, asPercentage=FALSE){
    standardGeneric("patternMatrix")
  }
)


setMethod("patternMatrix",
          signature(pattern = "DNAStringSet", windows = "DNAStringSet"),
          function(pattern, windows, cores = 1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in windows must have the same 
                    length!")
            }
            if(!(length(unique(width(pattern))) == 1)){
              stop("All sequences in pattern must have the same 
                    length!")
            }
            
            if(length(pattern)==1){
              # if there is only one pattern then create ScoreMatrix
	            mat <- .patternDNAStringSet_windowsDNAStringSet2matrix(pattern, windows, cores)
              print(as.data.frame(mat))
              print(new("ScoreMatrix",mat))
              return(new("ScoreMatrix",mat))
              
            }else{
              # if more than one pattern
              lmat <- lapply(1:length(pattern), 
                       function(i) .patternDNAStringSet_windowsDNAStringSet2matrix(pattern[i], 
                                                                                  windows, 
                                                                                  cores))
              return(new("ScoreMatrixList", lmat))
            }
})


setMethod("patternMatrix",
          signature(pattern = "DNAStringSet", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, min.score = NULL, cores=1, asPercentage=FALSE){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            if(!(length(unique(width(pattern))) == 1)){
              stop("All sequences in pattern must have the same 
                    length!")
            }
            
            #Error in .starfreeStrand(strand(names)) : 
            #cannot mix "*" with other strand values
            strand(windows) <- c("*")
            
            # Questions
            # ?trim windows that are out of chromosomes based on seqinfo from the BSgenome object
            seqinfo(windows) <- seqinfo(genome)[seqnames(seqinfo(windows))]
  
            #TODO: check if windows are not out of chromosomes?
            windows = getSeq(genome, windows)
            
            if(length(pattern)==1){              
              mat <- .patternDNAStringSet_windowsDNAStringSet2matrix(pattern, windows, cores)
              return(new("ScoreMatrix", mat))
            }else{
              # if more than one pattern
              lmat <- lapply(1:length(pattern), 
                             function(i) .patternDNAStringSet_windowsDNAStringSet2matrix(pattern, windows, cores))
              return(new("ScoreMatrixList", lmat))
            }
})


setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "DNAStringSet"),
          function(pattern, windows, min.score = NULL, asPercentage=FALSE){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                    length!")
            }
            if(class(pattern)=="list" &
               !all(sapply(1:length(patternlist), 
                          function(x) class(patternlist[[x]])=="matrix"))
              ){
              stop("Pattern argument contain not only matrices!")
            }
            

            if(length(pattern)==1 & class(patternlist[[1]]) == matrix){              
              mat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                            min.score, asPercentage)
              new("ScoreMatrix",mat)
            }else{
              # if more than one pattern
              lmat <- lapply(1:length(pattern), 
                             function(i) .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                                               min.score, 
                                                                               asPercentage))
              return(new("ScoreMatrixList", lmat))
            }
})



setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, 
                   min.score = NULL, asPercentage=FALSE){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            if(class(pattern)=="list" &
                 !all(sapply(1:length(patternlist), 
                             function(x) class(patternlist[[x]])=="matrix"))
            ){
              stop("Pattern argument contain not only matrices!")
            }
            
            #Error in .starfreeStrand(strand(names)) : 
            #cannot mix "*" with other strand values
            strand(windows) <- c("*")
            
            # Questions
            # ?trim windows that are out of chromosomes based on seqinfo from the BSgenome object
            seqinfo(windows) <- seqinfo(genome)[seqnames(seqinfo(windows))]
            # after shows info 
            #In valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE) :
	    #GRanges object contains 1 out-of-bound range located on sequence chr3.
      #Note that only ranges located on a non-circular sequence whose length
  #is not NA can be considered out-of-bound (use seqlengths() and
  #isCircular() to get the lengths and circularity flags of the underlying
  #sequences). You can use trim() to trim these ranges. See
  #?`trim,GenomicRanges-method` for more information.
            #windows = trim(windows)
  
            windows = getSeq(genome, windows)
  
            if(class(pattern)=="matrix"){              
                mat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                              min.score, asPercentage)
                new("ScoreMatrix",mat)
            }else{
                 print("DLACZEGOOO")
                 # if more than one pattern
                 lmat <- lapply(1:length(pattern), 
                                function(i) .patternPWM_windowsDNAStringSet2matrix(pattern[i], windows, 
                                                                                  min.score, asPercentage))
                 return(new("ScoreMatrixList", lmat))
            }
})





