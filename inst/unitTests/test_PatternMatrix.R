# ---------------------------------------------------------------------------- #
# tests for PatternMatrix function


test_PatternMatrix_matrix_DNAStringSet = function()
{
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

test_PatternMatrix_DNAStringSet_DNAStringSet = function()
{
  pattern <- DNAStringSet(c("AAA","CCC"))
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
  
}


#only check if it is somehow working
test_PatternMatrix_matrix_GRanges_BSgenome = function()
{

  #TODO: small BSgenome object?
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("BSgenome.Hsapiens.UCSC.hg19")
  library(BSgenome.Hsapiens.UCSC.hg19)
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
  #windows <- GRanges(seqnames = "chr3",
  #		ranges = IRanges(c(1,999999999), width = 10),
  #		strand = "+")
  
  # when windows contain only "N"
  windows <- GRanges(seqnames = "chr3",
                     ranges = IRanges(9001:9010, width = 10,
                                      names = head(letters,10)),
                     strand = "+")
  
  
}


test_PatternMatrix_DNAStringSet_GRanges_BSgenome = function()
{

  #TODO: small BSgenome object?
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("BSgenome.Hsapiens.UCSC.hg19")
  #library(BSgenome.Hsapiens.UCSC.hg19)
  #genome = BSgenome.Hsapiens.UCSC.hg19
  
  file <- system.file("extdata", "ce2chrM.fa.gz", package="BSgenome")
  forgeBSgenomeDataPkg(file)
  

  windows <- GRanges(seqnames = "chr3",
		ranges = IRanges(9001:90010, width = 10,
				 names = head(letters,10)),
		strand = "+")
  
  pattern <- DNAStringSet(c("AAA", "AAA", "AAA"))

  p = patternMatrix(pattern=pattern, windows=windows, genome=genome)
  
  # silent these warnings?
  #10: In .Call2("PWM_score_starting_at", pwm, subject, starting.at, base_codes,  :
  #'subject' contains letters not in [ACGT] ==> assigned weight 0 to them

  # what if coordinates out of chromosomes:
  #windows <- GRanges(seqnames = "chr3",
  #		ranges = IRanges(c(1,999999999), width = 10),
  #		strand = "+")
  
}


  
  
  
  
  
  
  
  
  
  
  






