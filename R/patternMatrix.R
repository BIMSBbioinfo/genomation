
#######################################
# S3 functions
#######################################

# It finds positions of sequence motif hits above
# a specified threshold in a list of sequences of the same length.
# Adapted from seqPattern::motifScanHits - 
# it returns vector of 0s and 1s instead of start positions
# of pattern occurences
# Args:
# pwm: matrix, seq: character
.scan.sequence.with.pwm <- function(pwm, seq, minScore){
  
  if(is.null(minScore)){
    # calculate scores themselves
    
    pwm.match <- matchPWM(pwm = pwm.corr, subject = seq, with.score=TRUE)
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
  
    pwm.match <- matchPWM(pwm = pwm.corr, subject = seq, min.score = score.threshold)
  
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

# Find positions of specified PWM 
# and calculate score matrix
# Args:
# pwm: matrix, windows: character
.patternPWM_windowsDNAStringSet2matrix = function(pattern, windows, 
                                                  min.score, asPercentage){
  print(".patternPWM_windowsDNAStringSet2matrix")
  print("min.score")
  print(min.score)
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


# Find positions of specified sequence patterns 
# and calculate score matrices (list of matrices)
# Args:
# pattern: character, windows: DNAStringSet
.patterncharacter_windowsDNAStringSet2matrix = function(pattern, windows){
  
  # getPatternOccurrenceList has nrCores and useMulticore args for paralelization,
  # it uses mclapply for every pattern
  # and before mclappling there is only one such line:
  # regionsSeq <- DNAStringSet(gsub("N", "+", regionsSeq)), 
  # where regionsSeq=windows 
  pat = getPatternOccurrenceList(windows, pattern)
  
  positions_of_01 = function(pat, pattern, windows){
    # gets the indeces where hits occur
    indeces=unlist(mapply(function(x,y) x:y, pat$position, pat$position+width(pattern[1])-1, SIMPLIFY=FALSE))
  
    # create output vector
    res=rep(0,length(windows[[1]]))
  
    # count how many times a PWM occurs in a given index
    ind.rle=rle(indeces)
  
    # get the indeces and assign the number of occurrence for that index
    res[ind.rle$values]=ind.rle$lengths
  
    return(res)
  }
  
  pat = pat[[1]]
  
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
#' Get scores that correspond to k-mer or PWM matrix occurence for bases in each window
#'
#'
#' @param pattern matrix (a PWM matrix) or a character vector.
#'                A matrix needs to have one row for each nucleotide ("A","C","G" and "T").
#'                IUPAC ambiguity codes can be used and will match any letter in the subject that is associated with the code.
#'                (pattern can not contain "N").
#' @param windows \code{\link{GRanges}} object or \code{\link{DNAStringSet}} object 
#'                that have the same width of ranges or sequences.
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
               asPercentage=FALSE,cores=1){
    standardGeneric("patternMatrix")
  }
)

setMethod("patternMatrix",
          signature(pattern = "character", windows = "DNAStringSet"),
          function(pattern, windows, cores = 1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in windows must have the same 
                   length!")
            }
            if(length(pattern)==1 & cores>1){
              stop("There is only one pattern provided. Setting number of cores to 1..")
              cores=1
            }

           if(length(nchar(pattern))==1){
                # if there is only one pattern then create ScoreMatrix
                mat <- .patterncharacter_windowsDNAStringSet2matrix(pattern, windows)
                
                return(new("ScoreMatrix",mat))
              
            }else{
               if(cores==1){
                 lmat <- lapply(1:length(pattern), 
                               function(i) .patterncharacter_windowsDNAStringSet2matrix(pattern[i],
                                                                                           windows))
               }else{
                 
                  lmat <- mclapply(pattern, 
                                 patterncharacter_windowsDNAStringSet2matrix,
                                 windows,
                                 mc.cores=cores)
               }
               
               return(new("ScoreMatrixList", lmat))
            }
})


setMethod("patternMatrix",
          signature(pattern = "character", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, min.score = NULL, asPercentage=FALSE, cores=1){
            
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
            
            #TODO: check if windows are not out of chromosomes?
            windows = getSeq(genome, windows)
            
            # call patternMatrix function
            # pattern: character, windows: DNAStringSet
            patternMatrix(pattern=pattern, windows=windows,
                          min.score=min.score, asPercentage=asPercentage,
                          cores=cores)
})


setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "DNAStringSet"),
          function(pattern, windows, min.score = NULL, asPercentage=FALSE, cores=1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            
            print("min.score")
            print(min.score)
            mat <- .patternPWM_windowsDNAStringSet2matrix(pattern, windows, 
                                                          min.score, asPercentage)
            return(new("ScoreMatrix",mat))
})


setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "GRanges", genome="BSgenome"),
          function(pattern, windows, genome, 
                   min.score = NULL, asPercentage=FALSE, cores=1){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
                   length!")
            }
            
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
            patternMatrix(pattern=pattern, windows=windows, 
                          min.score=min.score, asPercentage=asPercentage,
                          cores=1)
})   


#Error: in method for ‘patternMatrix’ with signature ‘pattern="list",windows="DNAStringSet"’: formal arguments in method and generic do not appear in the same order
# setMethod("patternMatrix",
#           signature(pattern = "list", windows = "DNAStringSet"),
#           function(pattern, windows, min.score = NULL, 
#                    asPercentage=FALSE, cores=1){
#             
#             if(class(pattern[[1]])=="matrix"){
#               
#               if(cores==1){
#                 lmat <- lapply(1:length(pattern), 
#                                function(i) patternMatrix(pattern=pattern[i], 
#                                                          windows=windows, 
#                                                          min.score=min.score, 
#                                                          asPercentage=asPercentage))
#               }else{
#                 lmat <- mclapply(pattern, 
#                                  patternMatrix,
#                                  windows=windows,
#                                  min.score=min.score,
#                                  asPercentage=asPercentage,
#                                  mc.cores=cores)
#               }
#               return(new("ScoreMatrixList",lmat))
#             }
#             
#             if(class(pattern[[1]])=="character"){
#               pattern = unlist(pattern)
#               patternMatrix(pattern, windows, min.score)
#             }
# })

#Error: in method for ‘patternMatrix’ with signature ‘pattern="list",windows="GRanges",genome="BSgenome"’: formal arguments in method and generic do not appear in the same order
# setMethod("patternMatrix",
#           signature(pattern = "list", windows = "GRanges", genome="BSgenome"),
#           function(pattern, windows, genome, min.score = NULL, 
#                    asPercentage=FALSE, cores=1){
#             
#             if(class(pattern[[1]])=="matrix"){
#               
#               if(cores==1){
#                 lmat <- lapply(1:length(pattern), 
#                                function(i) patternMatrix(pattern=pattern[i], 
#                                                          windows=windows, 
#                                                          genome=genome,
#                                                          min.score=min.score, 
#                                                          asPercentage=asPercentage))
#               }else{
#                 lmat <- mclapply(pattern, 
#                                  patternMatrix,
#                                  windows=windows,
#                                  genome=genome,
#                                  min.score=min.score,
#                                  asPercentage=asPercentage,
#                                  mc.cores=cores)
#               }
#               return(new("ScoreMatrixList",lmat))
#             }
#             
#             if(class(pattern[[1]])=="character"){
#               pattern = unlist(pattern)
#               patternMatrix(pattern, windows, genome=genome, 
#                             min.score=min.score, asPercentage=asPercentage)
#             }
# })


