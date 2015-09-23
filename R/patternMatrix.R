
# Notes:
# windows can have "N" but then score is set to 0 by
# patterns cannot contain "N"

# BSgenome:matchPattern

# Questions:
# if windows contain only N then return matrix with 0s or nothing (right now nothing)


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
              # they use parallel::mclapply for parallelization
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
              # calculate coverage of occurence of each pattern, 1-present, 0-absent
              ma <- lapply(1:length(windows), function(x) as.vector(coverage(reduce(pwm.match[[x]]))) )
              
              # convert list to matrix
              mat <- matrix(unlist(ma), ncol = length(ma[[1]]), byrow = TRUE)
              
            }
            
            return(mat)
}


#######################################
# S4 functions
#######################################


setGeneric(
  name="patternMatrix",
  def=function(pattern, windows, genome=NULL, min.score = NULL, cores=1, asPercentage=FALSE){
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
            
	      ma <- .patternDNAStringSet_windowsDNAStringSet2matrix(pattern, windows, cores)
              return(new("ScoreMatrix",mat))
              
            }else{
              #TODO
              print("ScoreMatrixList TODO")
            }
          }
)



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
            
            ma <- .patternDNAStringSet_windowsDNAStringSet2matrix(pattern, windows, cores)
            return(new("ScoreMatrix",mat))

            }
)


setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "DNAStringSet"),
          function(pattern, windows, min.score = NULL, asPercentage=FALSE){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
            length!")
            }
                        
            mat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                          min.score, asPercentage)
            new("ScoreMatrix",mat)

          }
)



setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, min.score = NULL, asPercentage=FALSE){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
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
            
            mat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                          min.score, asPercentage)
            new("ScoreMatrix",mat)
            
            }
)





