
#######################################
# S3 functions
#######################################

# It finds positions of sequence motif hits above
# a specified threshold in a list of sequences of the same length.
# Adapted from seqPattern::motifScanHits - 
# it returns vector of 0s and 1s instead of start positions
# of pattern occurrences
# Args:
# pwm: matrix, seq: character
.scan.sequence.with.pwm <- function(seq, pwm, minScore){
  
  if(is.null(minScore)){
    # calculate scores themselves
    
    # suppressWarnings is done because matchPWM function gives warnings 
    # if you there is "N" nucleotide within sequence (can be only [A,T,C,G])
    # It assigns weight 0 to them.
    pwm.match <- suppressWarnings(matchPWM(pwm = pwm, subject = seq, with.score=TRUE))
    return(mcols(pwm.match)$score)
    
  }else{
    # calculate 0 and 1s
    
    # if minScore is a character then convert it to a number 
    # (e.g. 80% to 0.8)
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
  
    # suppressWarnings is done because matchPWM function gives warnings 
    # if you there is "N" nucleotide within sequence (can be only [A,T,C,G]).
    # It assigns weight 0 to them.
    pwm.match <- suppressWarnings(matchPWM(pwm = pwm.corr, subject = seq, min.score = score.threshold))
  
    # if PWM does not match any position in the sequence then return only 0s
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
}

# Find positions of specified PWM and calculate score matrix
# Args:
# pwm: matrix, windows: character
.patternPWM_windowsDNAStringSet2matrix = function(pattern, windows, 
                                                  min.score, asPercentage,
                                                  cores){
  if(is.null(min.score)){
    # for each window calculate vector of scores themselves
    
    # suppressWarnings is done because motifScanScores function gives warnings 
    # if you there is "N" nucleotide within sequence (can be only [A,T,C,G])
    # It assigns weight 0 to them.
    # motifScanScores is fast
    if(cores==1){
      mat <- suppressWarnings(motifScanScores(regionsSeq = windows,
                                             motifPWM = pattern,
                                             asPercentage = asPercentage))
    }else{
      mat <- mclapply(windows, 
                      motifScanScores,
                      motifPWM = pattern,
                      asPercentage = asPercentage,
                      mc.cores=cores)
      mat <- matrix(unlist(mat), ncol = length(mat[[1]]), byrow = TRUE)
    }

  # min.score!=NULL
  }else{
    pwm.match <- mclapply(windows, 
                          .scan.sequence.with.pwm,
                          pwm=pattern,
                          minScore=min.score,
                          mc.cores=cores
                          )
    # convert list to matrix
    mat <- matrix(unlist(pwm.match), ncol = length(pwm.match[[1]]), byrow = TRUE)
  }
  return(mat)
}

.getPatternOccurrenceList.wrapper = function(i, windows, pattern){
  # [i] instead of [[i]]
  pat <- getPatternOccurrenceList(windows[i], pattern)
  pat <- pat[[1]]
  pat$sequence <- rep(i, nrow(pat))
  pat
}

# Find positions of specified sequence patterns 
# and calculate score matrices (list of matrices)
# Args:
# pattern: character, windows: DNAStringSet
.patterncharacter_windowsDNAStringSet2matrix = function(pattern, windows, cores){

  # getPatternOccurrenceList has nrCores and useMulticore args for paralelization,
  # it uses mclapply for every pattern
  # and before mclappling there is only one such line:
  # regionsSeq <- DNAStringSet(gsub("N", "+", regionsSeq)), 
  # where regionsSeq=windows 
  if(cores==1){
    pat = getPatternOccurrenceList(windows, pattern)
    pat = pat[[1]]
  }else{
    pat = mclapply(1:length(windows),
                   .getPatternOccurrenceList.wrapper,
                   windows=windows,
                   pattern=pattern,
                   mc.cores=cores)
    
    # data.table:rbindlist
    # Same as do.call("rbind", l) on data.frames, but much faster
    pat = rbindlist(pat)
  }
  
  positions_of_01 = function(pat, pattern, windows){
    # gets the indeces where hits occur
    indeces=unlist(mapply(function(x,y) x:y, 
                          pat$position, pat$position+width(pattern[1])-1, SIMPLIFY=FALSE))
  
    # create output vector
    res=rep(0,length(windows[[1]]))
  
    # count how many times a PWM occurs in a given index
    ind.rle=rle(indeces)
  
    # get the indeces and assign the number of occurrence for that index
    res[ind.rle$values]=ind.rle$lengths
  
    return(res)
  }

  # number windows from result of getPatternOccurrenceList
  pat.windows = unique(pat$sequence)
  
  # for each window calculate vector of 1 and 0s
  pwm.match <- lapply(1:length(windows), function(i){
    if(i %in% pat.windows){
      positions_of_01(pat[pat$sequence==i,], pattern, windows)
    }else{
      # if pattern did not match to window
      rep(0,length(windows[[1]]))
    }
  })
  # convert list to matrix
  mat <- matrix(unlist(pwm.match), ncol = length(pwm.match[[1]]), byrow = TRUE)
}


#######################################
# S4 functions
#######################################

# ---------------------------------------------------------------------------- #
#' Get scores that correspond to k-mer or PWM matrix occurrence for bases in each window
#'
#' The function produces a base-pair resolution matrix or matrices of scores that correspond to 
#' k-mer or PWM matrix occurrence over predefined windows that  have equal width.
#' It finds either positions of pattern hits above a specified threshold and 
#' creates score matrix filled with 1 (presence of pattern) and 0 (its absence) or
#' matrix with scores themselves.
#' If pattern is a character of length 1 or PWM matrix then the function returns 
#' a ScoreMatrix object, if character of length more than 1 or list of PWMs 
#' then ScoreMatrixList.
#'
#' @param pattern matrix (a PWM matrix), list of matrices or a character vector of length 1 or more.
#'                A matrix is a PWM matrix that needs to have one row for each nucleotide 
#'                ("A","C","G" and "T" respectively).
#'                IUPAC ambiguity codes can be used and it will match any letter 
#'                in the subject that is associated with the code.
#' @param windows \code{\link{GRanges}} object or \code{\link{DNAStringSet}} object 
#'                that have equal width of ranges or sequences.
#' @param genome \code{\link{BSgenome}} object
#' @param min.score numeric or character indicating minimum score to count a match. 
#'                  It can be given as a character string containing
#'                  a percentage of the highest possible score or a single number
#'                  (by default "80\%" or 0.8). If min.score is set to NULL
#'                  then \code{patternMatrix} returns scores themselves (default).
#' @param asPercentage boolean telling whether scores represent percentage of the maximal 
#'                     motif PWM score (default: TRUE) or raw scores (FALSE).             
#' @param cores the number of cores to use (default: 1). It is supported only on Unix-like platforms.
#'              
#' @details
#' \code{patternMatrix} is based on functions from the seqPattern package:
#' getPatternOccurrenceList function to find position of pattern that is a character vector in 
#' a list of sequences (a DNAStringSet object)
#' and adapted function motifScanHits to find pattern that is a PWM matrix
#' in sequences (a DNAStringSet object).
#' 
#' If cores > 1 is provided then for every window occurrence of pattern is counted in paralallel.
#' 
#' @examples
#' library(Biostrings)
#' 
#' # consensus sequence of the ctcf motif
#' motif = "CCGCGNGGNGGCAG"
#' # Creates 10 random DNA sequences
#' seqs = sapply(1:10, 
#'        function(x) paste(sample(c("A","T","G","C"), 180, replace=TRUE), collapse=""))
#' windows = DNAStringSet(seqs)
#' p = patternMatrix(pattern=motif, windows=windows, min.score=0.8)
#' p
#' 
#' @return returns a \code{scoreMatrix} object or a \code{scoreMatrixList} object
#' 
#' @seealso \code{\link{ScoreMatrix}}, \code{\link{ScoreMatrixList}}
#' @docType methods
#' @rdname patternMatrix-methods           
#' @export
setGeneric(
  name="patternMatrix",
  def=function(pattern, windows, 
               genome=NULL, min.score = 0.8, 
               asPercentage=FALSE, cores=1){
    standardGeneric("patternMatrix")
  }
)

#' @aliases patternMatrix,character,DNAStringSet-method
#' @rdname patternMatrix-methods
#' @usage  \\S4method{patternMatrix}{character,DNAStringSet}(pattern, windows,
#'                                                           asPercentage, cores)
setMethod("patternMatrix",
          signature(pattern = "character", windows = "DNAStringSet"),
          function(pattern, windows, cores = 1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in windows must have the same 
                   length!")
            }
            if(length(pattern)==1 & cores>1){
              warning("There is only one pattern provided. Setting number of cores to 1..")
              cores=1
            }

           if(length(nchar(pattern))==1){
                # if there is only one pattern then create ScoreMatrix
                mat <- .patterncharacter_windowsDNAStringSet2matrix(pattern, windows, cores=cores)
                return(new("ScoreMatrix",mat))
              
            }else{
                lmat <- lapply(1:length(pattern), 
                             function(i) .patterncharacter_windowsDNAStringSet2matrix(pattern[i],
                                                                                      windows,
                                                                                      cores=cores))
                return(new("ScoreMatrixList", lmat))
            }
})

#' @aliases patternMatrix,character,GRanges,BSgenome-method
#' @rdname patternMatrix-methods
#' @usage  \\S4method{patternMatrix}{character,GRanges,BSgenome}(pattern, windows, genome,
#'                                                               cores)
setMethod("patternMatrix",
          signature(pattern = "character", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, cores=1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            
            #Error in .starfreeStrand(strand(names)) : 
            #cannot mix "*" with other strand values
            strand(windows) <- c("*")
            
            seqinfo(windows) <- seqinfo(genome)[seqnames(seqinfo(windows))]
            
            windows = getSeq(genome, windows)
            
            # call patternMatrix function
            # pattern: character, windows: DNAStringSet
            patternMatrix(pattern=pattern, windows=windows,
                          cores=cores)
})

#' @aliases patternMatrix,matrix,DNAStringSet-method
#' @rdname patternMatrix-methods
#' @usage  \\S4method{patternMatrix}{matrix,DNAStringSet}(pattern, windows,
#'                                                        min.score, asPercentage, 
#'                                                        cores)
setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "DNAStringSet"),
          function(pattern, windows, min.score = 0.8, asPercentage=FALSE, cores=1){
                        
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            
            mat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                          min.score, asPercentage,
                                                          cores=cores)
            return(new("ScoreMatrix",mat))
})

#' @aliases patternMatrix,matrix,GRanges,BSgenome-method
#' @rdname patternMatrix-methods
#' @usage  \\S4method{patternMatrix}{matrix,GRanges,BSgenome}(pattern, windows, genome,
#'                                                            min.score, asPercentage, 
#'                                                            cores)
setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, 
                   min.score = 0.8, asPercentage=FALSE, cores=1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            
            strand(windows) <- c("*")
            
            seqinfo(windows) <- seqinfo(genome)[seqnames(seqinfo(windows))]
            # after it can show info 
            #....
            #Note that only ranges located on a non-circular sequence whose length
            #is not NA can be considered out-of-bound (use seqlengths() and
            #isCircular() to get the lengths and circularity flags of the underlying
            #sequences). You can use trim() to trim these ranges. See
            #?`trim,GenomicRanges-method` for more information.
            #windows = trim(windows)
            
            windows = getSeq(genome, windows)
            
            # call patternMatrix function
            # pattern: matrix, windows: DNAStringSet
            patternMatrix(pattern=pattern, windows=windows, 
                          min.score=min.score, asPercentage=asPercentage,
                          cores=1)
})   

#' @aliases patternMatrix,list,DNAStringSet-method
#' @rdname patternMatrix-methods
#' @usage  \\S4method{patternMatrix}{list,DNAStringSet}(pattern, windows,
#'                                                      min.score, asPercentage, 
#'                                                      cores)
setMethod("patternMatrix",
          signature(pattern = "list", windows = "DNAStringSet"),
          function(pattern, windows, min.score = 0.8, 
                   asPercentage=FALSE, cores=1){
            
            if(inherits(pattern[[1]], "matrix")){
              lmat <- lapply(1:length(pattern), 
                            function(i) patternMatrix(pattern=pattern[[i]], 
                                                      windows=windows, 
                                                      min.score=min.score, 
                                                      asPercentage=asPercentage,
                                                      cores=cores))
              return(new("ScoreMatrixList",lmat))
            }
            
            if(is.character(pattern[[1]])){
              pattern = unlist(pattern)
              patternMatrix(pattern, windows, min.score)
            }
})

#' @aliases patternMatrix,list,GRanges,BSgenome-method
#' @rdname patternMatrix-methods
#' @usage  \\S4method{patternMatrix}{list,GRanges,BSgenome}(pattern, windows, genome, 
#'                                                          min.score, asPercentage, 
#'                                                          cores)
setMethod("patternMatrix",
          signature(pattern = "list", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, min.score = 0.8, 
                   asPercentage=FALSE, cores=1){
                   
            if(inherits(pattern[[1]],"matrix")){
              lmat <- lapply(1:length(pattern), 
                            function(i) patternMatrix(pattern=pattern[i], 
                                                      windows=windows, 
                                                      genome=genome,
                                                      min.score=min.score, 
                                                      asPercentage=asPercentage,
                                                      cores=cores))
              return(new("ScoreMatrixList",lmat))
            }
            
            if(is.character(pattern[[1]])){
              pattern = unlist(pattern)
              patternMatrix(pattern, windows, genome=genome, 
                            min.score=min.score, asPercentage=asPercentage)
            }
})


