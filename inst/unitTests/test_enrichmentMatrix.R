# ---------------------------------------------------------------------------- #
# enrichmentMatrix()
# ---------------------------------------------------------------------------- #
bam.file <- system.file('unitTests/test.bam', package='genomation')
windows <- GRanges(rep(c(1,2), each=2), IRanges(rep(c(1, 2), times=2), width=5))
IP <- ScoreMatrix(target=bam.file,windows=windows, type='bam')

bam.file <- system.file('unitTests/test.bam', package='genomation')
control <- ScoreMatrix(target=bam.file, windows=windows, type='bam') * 0.002

# IP - ScoreMatrix, control - ScoreMatrix
esm_1 <- enrichmentMatrix(IP, control)
esm_2 <-log2((IP + 1) / (control + 1))
checkEquals(esm_1, esm_2)

# IP - ScoreMatrixList, control - ScoreMatrixList
sml_IP <- ScoreMatrixList(list(IP1 = IP, IP2 = IP))
sml_control <- ScoreMatrixList(list(c1 = control, c2 = control))
esml_1 <- enrichmentMatrix(sml_IP, sml_control)
es1 <- log2((IP + 1) / (control + 1))
esml_2 <- ScoreMatrixList(list(IP1 = es1, IP2 = es1))
checkEquals(esml_1, esml_2)

# IP - ScoreMatrixList, control - ScoreMatrix
sml_IP <- ScoreMatrixList(list(IP1 = IP, IP2 = IP))
esml_3 <- enrichmentMatrix(sml_IP, control)
es2 <- log2((IP + 1) / (control + 1))
esml_4 <- ScoreMatrixList(list(IP1 = es1, IP2 = es1))
checkEquals(esml_3, esml_4)

# error for unequal list lengths
checkException(enrichmentMatrix(ScoreMatrixList(list(IP, IP, IP)), sml_control), silent=TRUE)
