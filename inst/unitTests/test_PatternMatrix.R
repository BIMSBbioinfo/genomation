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


test_PatternMatrix_matrix_BSgenome = function()
{
  
  #TODO
  
  source("https://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  genome = BSgenome.Hsapiens.UCSC.hg19
  
  seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
  gr <-
    GRanges(seqnames =
              Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges = IRanges(
              1:10, width = 10:1, names = head(letters,10)),
            strand = Rle(
              strand(c("-", "+", "*", "+", "-")),
              c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10),
            seqinfo=seqinfo(genome))
  
  seqinfo(gr) <- seqinfo(genome)
  
  #Error in .starfreeStrand(strand(names)) : 
  #cannot mix "*" with other strand values
  strand(gr) <- c("*")
  
  
  
}