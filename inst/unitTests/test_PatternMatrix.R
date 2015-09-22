# ---------------------------------------------------------------------------- #
# test for PatternMatrix function

#only check if is somehow working

test_PatternMatrix_matrix_DNAStringSet = function()
{
  patterns <- DNAStringSet(c("GCT", "GGT", "GCA"))
  pwm <- PWM(patterns)
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),DNAString("AAAGCTAAAGGTAAAGCAAAA") ,
                    DNAString("AAAGCTAAAGGTAAAGCAAAA") ), "DNAStringSet")
  p = patternMatrix(pattern=pwm, windows=windows)
  
  
}


test_PatternMatrix_DNAStringSet_DNAStringSet = function()
{
  pattern <- DNAStringSet(c("GCT"))
  windows = as(list(DNAString("AAAGCTAAAGGTAAAGCAAAA"),DNAString("AAAGCTAAAGGTAAAGCAAAA") ,
                    DNAString("AAAGCTAAAGGTAAAGCAAAA") ), "DNAStringSet")
  p = patternMatrix(pattern=pattern, windows=windows)
  
}