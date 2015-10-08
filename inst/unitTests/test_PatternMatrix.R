# ---------------------------------------------------------------------------- #
# tests for PatternMatrix function


test_PatternMatrix_matrix_DNAStringSet = function()
{
  # one pwm matrix
  patterns <- DNAStringSet(c("AAA", "AAA", "AAA"))
  pwm <- PWM(patterns)
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),
		    DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA")), "DNAStringSet")
  p1 = patternMatrix(pattern=pwm, windows=windows, min.score=1)
  mat1=as((as.matrix(windows)=="A")*1, "ScoreMatrix")
  
  checkEquals(p1, mat1)
}

hits <- matchPWM(pwm, windows[[1]], with.score=TRUE)

test_PatternMatrix_character_DNAStringSet = function()
{
  # one pattern
  pattern <- "AAA"
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),
		    DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA")), "DNAStringSet")
  p2 = patternMatrix(pattern=pattern, windows=windows)
  mat1=as((as.matrix(windows)=="A")*1, "ScoreMatrix")
  
  checkEquals(p2, mat1) 

  # Questions
  #pattern <- DNAStringSet(c("AAA"))
  #windows = as(list(DNAString("NNNNNNNNNNNNNN"),
  #		    DNAString("NNNNNNNNNNNNNN")), "DNAStringSet")
  #p2 = patternMatrix(pattern=pattern, windows=windows)
  #mat1=as((as.matrix(windows)=="A")*1, "ScoreMatrix")
  #checkEquals(p2, mat1) 
  
  
  # more than one pattern
  pattern <- c("AAA","AAA")
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA")), "DNAStringSet")
  p3 = patternMatrix(pattern=pattern, windows=windows)
  mat3=as(list((as.matrix(windows)=="A")*1, (as.matrix(windows)=="A")*1), 
          "ScoreMatrixList")
  checkEquals(p3, mat3) 
  
}

#TODO: how to create small BSgenome object?
test_PatternMatrix_matrix_GRanges_BSgenome = function()
{
  
  # one PWM matrix

  #TODO: how to get small BSgenome object?
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("BSgenome.Hsapiens.UCSC.hg19")
  #library(BSgenome.Hsapiens.UCSC.hg19)
  genome = BSgenome.Hsapiens.UCSC.hg19

  windows <- GRanges(seqnames = "chr3",
		ranges = IRanges(9001:9010, width = 10,
				 names = head(letters,10)),
		strand = "+")
  
  patterns <- DNAStringSet(c("AAA", "AAA", "AAA"))
  pattern <- PWM(patterns)

  p = patternMatrix(pattern=pattern, windows=windows, genome=genome)
  
  # silent these warnings?
  #10: In .Call2("PWM_score_starting_at", pwm, subject, starting.at, base_codes,  :
  #'subject' contains letters not in [ACGT] ==> assigned weight 0 to them

  # what if coordinates out of chromosomes:
  # windows <- GRanges(seqnames = "chr3",
  #		ranges = IRanges(c(1,999999999), width = 10),
  #		strand = "+")
  
  # when windows contain only "N"
  #windows <- GRanges(seqnames = "chr3",
  #                   ranges = IRanges(9001:9010, width = 10,
  #                                    names = head(letters,10)),
  #                   strand = "+")
}


test_PatternMatrix_character_GRanges_BSgenome = function()
{
  
  genome = BSgenome.Hsapiens.UCSC.hg19
  
  windows <- GRanges(seqnames = "chr3",
		ranges = IRanges(9001:90010, width = 10,
				 names = head(letters,10)),
		strand = "+")
  
  # one pattern
  pattern <- "AAA"
  p = patternMatrix(pattern=pattern, windows=windows, genome=genome)
  
  # more than one pattern
  pattern <- c("AAA","AAA")
  p = patternMatrix(pattern=pattern, windows=windows, genome=genome)
  
  # silent these warnings?
  #10: In .Call2("PWM_score_starting_at", pwm, subject, starting.at, base_codes,  :
  #'subject' contains letters not in [ACGT] ==> assigned weight 0 to them

  # what if coordinates out of chromosomes:
  # windows <- GRanges(seqnames = "chr3",
  #		ranges = IRanges(c(1,999999999), width = 10),
  #		strand = "+")
}


  
  
  
  
  
  
  
  
  
  
  






