
# Notes:
# windows can have "N" but then score is set to 0 by
# patterns can not contain "N" ->  getPatternOccurrenceList returns NULL

# quantile(p@.Data, probs = c(0, 0.25, 0.5, 0.75, 0.9))

# BSgenome:matchPattern

# Questions:
# if windows contain only N then return matrix with 0s or nothing (right now nothing)

#For PWMs, we can either return occurence if the user provides a score cutoff, 
#or return scores themselves, I think the other function, motifScanHits() is capable of doing both
# TODO:  some basic benchmarks for 1000,5000,10000 windows

#pattern doesnt have to have the same length

#######################################
# S3 functions
#######################################

# It finds positions of sequence motif hits above
# a specified threshold in a list of sequences of the same length.
# Adapted from seqPattern::motifScanHits - 
# it returns granges instead of start(granges).
.scan.sequence.with.pwm <- function(pwm, seq, minScore){
  
  # convert minScore character to number (e.g. 80% to 0.8)
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
  
  # if pwm does not match to any position in the sequence
  # return only 0s
  if(length(pwm.match)==0){
    return(rep(0,length(seq)))
  }
  
  # gets the indeces where hits occur
  indeces=unlist(mapply(function(x,y) x:y, start(pwm.match),end(pwm.match), SIMPLIFY=FALSE))
  
  # create output vector
  res=rep(0,length(seq))
  
  # count how many times a PWM occurs in a given index
  ind.rle=rle(indeces)
  
  # get the indeces and assign the number of occurrence for that index
  res[ind.rle$values]=ind.rle$lengths
  
  return(res)
}


# Find positions of specified sequence patterns 
# and calculate score matrices (list of matrices)
.patternDNAStringSet_windowsDNAStringSet2list = function(pattern, windows, cores){
  
  if(class(pattern)=="list"){
    # because getPatternOccurrenceList takes a character vector as pattern arg
    pattern <- unlist(pattern) 
  }
  
  # Calculate start positions of specified sequence patterns
  if(length(nchar(pattern)) > 1){
    # create ScoreMatrixList if more than one pattern
    if(cores>1){
      # seqPattern::getPatternOccurrenceList uses parallel::mclapply for parallelization
      pat = getPatternOccurrenceList(windows, pattern,
                                     useMulticore=TRUE, nrCores=cores)
    }
  }else{
    # create ScoreMatrixList if >1 pattern or ScoreMatrix if only 1 pattern
    pat = getPatternOccurrenceList(windows, pattern)
  }  
  
  # pat is always a list       
  # check if patterns map into sequences that contain only "N" 'nucleotide'
  isnull <- sapply(1:length(pat), function(x) is.null(pat[x]))
  if(all(isnull)){
    stop(paste0("All patterns map into windows that contain letters not in [ACGT],",
                " e.g. they contain only 'N'"))
  }else if(any(isnull)){
    warning(paste0("All windows in ", paste(which(isnull), sep=","),
                   " element of the list contain letters not in [ACGT], e.g. they contain only 'N'.",
                   " Selecting only those that are valid")) #TODO
    pat <- pat[which(!isnull)]
  }
  
  calc.OccurenceCov = function(pat, pattern){
    # convert data.frame with coordinates of the positions within which sequence pattern occurs to GRanges
    # TODO: sprawdzic czy napewno dlugosc tego gr sie zgadza z dlugoscia patternu
    gr <- GRanges(seqnames = "chrN",
                  ranges = IRanges(start=pat$position, 
                                   width=width(pattern)),
                  strand = "*",
                  sequence = pat$sequence)
    #seqinfo = Seqinfo(seqnames="chrN",
    #                   seqlengths=max(end(windows)))) # po co to??
    
    # some of windows are not included, because pattern did not map into them
    seqq <- as.factor(gr$sequence)
    levels(seqq) <- as.character(1:length(windows))
    grl <- split(gr, seqq, drop=FALSE)
    
    # calculate coverage of occurence of each pattern per base pair, 1-present, 0-absent
    grl.reduce <- reduce(grl)
    ma <- lapply(1:length(grl.reduce), 
                 function(x) coverage(grl.reduce[x])[[1]] )
    # convert list to matrix
    mat <- matrix(unlist(ma), ncol = length(ma[[1]]), byrow = TRUE)
    return(mat)
  }
  
  # for each pattern calculate .. TODO
  pattern = list(pattern)
  #return(lapply(1:length(pat), function(i) calc.OccurenceCov(pat[[i]], pattern[[i]])))
}


.patternPWM_windowsDNAStringSet2matrix = function(pattern, windows, 
                                                  min.score, asPercentage){
  if(is.null(min.score)){
    # for each window calculate vector of scores themselves
    mat <- motifScanScores(regionsSeq = windows,
                           motifPWM = pattern,
                           asPercentage = asPercentage)
  }else{
    # for each window calculate vector of 1 and 0s
    pwm.match <- lapply(1:length(windows), function(x){
      .scan.sequence.with.pwm(pwm = pattern, seq = windows[[x]], minScore=min.score) 
    })
    
    # convert list to matrix
    mat <- matrix(unlist(pwm.match), ncol = length(pwm.match[[1]]), byrow = TRUE)
  }
  return(mat)
}



#######################################
# S4 functions
#######################################

# ---------------------------------------------------------------------------- #
#' Get scores that correspond to k-mer or PWM matrix occurence for bases in each window
#'
#' we are looking for occurrence, we aim to get a scoreMatrix with 1s and 0s for k-mers.
#' The function first look for and calculates
#' The function firsts bins each window to equal number of bins, and calculates
#' the a summary matrix for scores of each bin (currently, mean, max and min supported)
#' A patternMatrix object can be used to ..
#' \code{windows} can be a predefined region such as CpG islands or gene bodies that have the same width.
#'
#' @param pattern matrix (a PWM matrix), a DNAStringSet object or list of matrices or DNAStringSet objects. 
#'                A matrix has to have one row for each nucleotide ("A","C","G" and "T").
#'                IUPAC ambiguity codes can be used and will match any letter in the subject that is associated with the code.
#'                (can not contain "N").
#' @param windows \code{\link{GRanges}} object 
#' @param genome \code{\link{BSgenome}} object
#' @param min.score numeric or string indicating minimum score to count a match. 
#'                  It can be given as a character string containing
#'                  a percentage of the highest possible score or as a single number 
#'                  (as default "80\%"). If min.score is not provided or is set to NULL
#'                  then \code{patternMatrix} returns scores themselves.
#' @param asPercentage boolean telling whether scores represent percentage of the maximal 
#'                     motif PWM score (default: TRUE) or raw scores (FALSE).               
#' @param cores the number of cores to use (default: 1). It is supported only on Unix-like platforms.
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
            if(length(pattern)==1 & cores>1){
              stop("There is only one pattern provided. Setting number of cores to 1..")
              cores=1
            }
            
            #TODO
            if(cores==1){
              # if there is only one pattern then create ScoreMatrix
              mat <- .patternDNAStringSet_windowsDNAStringSet2matrix(pattern, windows, cores)
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
            
            # call patternMatrix function
            # pattern: DNAStringSet, windows: DNAStringSet
            pattern(pattern, windows, min.score)
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
             
            #returns matrix or list of matrices
            mat.or.lmat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                          min.score, asPercentage)
            
            if(class(mat.or.lmat)=="matrix"){
              return(new("ScoreMatrix",mat.or.lmat))
            }
            return(new("ScoreMatrixList",mat.or.lmat))
            
            })

setMethod("patternMatrix",
          signature(pattern = "list", windows = "DNAStringSet"),
          function(pattern, windows, min.score = NULL, asPercentage=FALSE){
            
            if(class(pattern[[1]])=="matrix"){
              
            }
            if(class(pattern[[1]])=="character" | class(pattern[[1]])=="DNAString"){
              pattern = DNAStringSet(unlist(pattern))
            }
            # call patternMatrix function
            sml = patternMatrix(pattern, windows, min.score)
            return(new("ScoreMatrixList",sml))
          })


setMethod("patternMatrix",
          signature(pattern = "list", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, 
                   min.score = NULL, asPercentage=FALSE,
                   cores=1){
            
            if(class(pattern[[1]]) %in% c("character", "DNAString", "matrix")){
              stop("Pattern has to be a list of character, DNAString or matrix objects")
            }
            
            # for each PWM calculate score matrix separately..
            if(class(pattern[[1]])=="matrix"){
              
              if(cores==1){
                sml <- lapply(1:length(pattern), 
                              function(i) patternMatrix(pattern=pattern[[i]], 
                                                        windows=windows, 
                                                        min.score=min.score,
                                                        asPercentage=asPercentage))
              }if else(cores > 1){
                sml <- mclapply(pattern,
                                patternMatrix,
                                  windows=windows,
                                  min.score=min.score,
                                  asPercentage=asPercentage,
                                  mc.cores = cores
                                )
              }else{
                stop("Numer of cores has to be >= 1.")
              }
              return(new("ScoreMatrixList",sml))
            }
            
            
            if(class(pattern[[1]])=="character"){
              pattern = DNAStringSet(unlist(pattern))
            }else if(class(pattern[[1]])=="DNAString"){
              pattern = DNAStringSet(pattern)
            }
            # call patternMatrix function.
            # if pattern is a DNAStringSet, then 
            # seqPattern::getPatternOccurrenceList function is used
            # and for parallelization it uses parallel::mclapply
            sml = pattern(pattern, windows, min.score)  
            return(new("ScoreMatrixList",sml))
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
            
            # call patternMatrix function
            # pattern: matrix, windows: DNAStringSet
            pattern(pattern, windows, min.score)
})

