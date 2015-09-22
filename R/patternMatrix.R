#######################################
# S3 functions
#######################################

# adated from seqPattern::motifScanHits
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


#######################################
# S4 functions
#######################################


#useMulticore  
#Logical, should multicore be used. useMulticore = TRUE is supported only on Unix-like platforms.
#nrCores	
#Number of cores to use when useMulticore = TRUE. Default value NULL uses all detected cores.
#it uses mclapply function

# i dont use getPatternOccurrenceList, i use it actually.
# because 
# https://github.com/Bioconductor-mirror/seqPattern/blob/master/R/MotifScanningMethods.R

#CHeck as percentage


setGeneric(
  name="patternMatrix",
  def=function(pattern, windows, min.score = NULL, cores=1, asPercentage=FALSE){
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
                        
            if(cores>1){
              # they use parallel::mclapply for parallelization
              pat = getPatternOccurrenceList(windows, pattern,
                                             useMulticore=TRUE, nrCores=cores)
            }else{
              pat = getPatternOccurrenceList(windows, pattern)
            }
            
            if(length(pattern)==1){
              
              pat = pat[[1]]
              gr <- GRanges(seqnames = "chrN",
                            ranges = IRanges(start=pat$position, width=length(pattern)),
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
              
              return(new("ScoreMatrix",mat))
            }else{
              #TODO
              print("ScoreMatrixList TODO")
            }
          }
)



setMethod("patternMatrix",
          signature(pattern = "matrix", windows = "DNAStringSet"),
          function(pattern, windows, min.score = NULL, asPercentage=FALSE){
            
            if(!(length(unique(width(windows))) == 1)){
              stop("All sequences in the input DNAStringSet must have the same 
            length!")
            }
                        
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
            
            new("ScoreMatrix",mat)

          }
)





