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
  
  #more than one pwm matrix
  patterns <- DNAStringSet(c("AAA", "AAA", "AAA"))
  pwmlist <- list(PWM(patterns),PWM(patterns))
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA")), "DNAStringSet")
  p1.1 = patternMatrix(pattern=pwmlist, windows=windows, min.score=1)
  sm1.1=as((as.matrix(windows)=="A")*1, "ScoreMatrix")
  sml1.1 <- as(list(sm1.1, sm1.1), "ScoreMatrixList") 
  checkEquals(p1.1, sml1.1)

}

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

  # when windows have only "N"
  pattern <- "AAA"
  windows = as(list(DNAString("NNNNNNNNNNNNNN"),
  		    DNAString("NNNNNNNNNNNNNN")), "DNAStringSet")
  p3 = patternMatrix(pattern=pattern, windows=windows)
  mat2=as((as.matrix(windows)=="A")*0, "ScoreMatrix")
  checkEquals(p3, mat2) 
  
  
  # more than one pattern
  pattern <- c("AAA","AAA")
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA"),
                    DNAString("AAAGCTAAAGGTAAAGCAAAA")), "DNAStringSet")
  p4 = patternMatrix(pattern=pattern, windows=windows)
  sml1=as(list((as.matrix(windows)=="A")*1, (as.matrix(windows)=="A")*1), 
          "ScoreMatrixList")
  checkEquals(p4, sml1) 
  
}
