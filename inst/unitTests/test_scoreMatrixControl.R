# ---------------------------------------------------------------------------- #
# scoreMatrixControl
# ---------------------------------------------------------------------------- #

bam.file <- system.file('unitTests/test.bam', package='genomation')
windows <- GRanges(rep(c(1,2), each=2), IRanges(rep(c(1, 2), times=2), width=5))
IP <- ScoreMatrix(target=bam.file,windows=windows, type='bam')

bam.file <- system.file('unitTests/test.bam', package='genomation')
control <- ScoreMatrix(target=bam.file, windows=windows, type='bam') * 0.002

# ------------------------------- #
# return ScoreMatrixControl
smc1 <- ScoreMatrixControl(IP, control)
smc2 <- as(IP, "ScoreMatrixControl")
smc2@control <- control
checkEquals(smc1, smc2)

# error for unequal matrix sizes
checkException(ScoreMatrixControl(IP, control[1:2,]), silent=TRUE)

# ---------------------------------------------------------------------------- #
# scoreMatrixListControl
# ---------------------------------------------------------------------------- #
sml_IP <- ScoreMatrixList(list(IP1 = IP, IP2 = IP))
sml_control <- ScoreMatrixList(list(c1 = control, c2 = control))

smlc1 <- ScoreMatrixListControl(sml_IP, sml_control)

## return ScoreMatrixListControl
# when control is a ScoreMatrixList object
smlc2 <- as(list(IP1 = smc1, IP2 = smc1), "ScoreMatrixListControl")
checkEquals(smlc1, smlc2)

## when control is a ScoreMatrix object
smlc3 <- ScoreMatrixListControl(sml_IP, control)
checkEquals(smlc3, smlc1)

# error for unequal list lengths
checkException(ScoreMatrixListControl(ScoreMatrixList(list(IP, IP, IP)), sml_control), silent=TRUE)

# ---------------------------------------------------------------------------- #
# c() function
# ---------------------------------------------------------------------------- #

# label names are given by a user
smlc_1 <- c(smlc1, s1 = smc1, s2 = smc1)
sml_IP2 <- ScoreMatrixList(list(IP1 = IP, IP2 = IP, s1 = IP, s2 = IP))
sml_control2 <- ScoreMatrixList(list(control, control, control, control))
smlc_2 <- ScoreMatrixListControl(sml_IP2, sml_control2)
checkEquals(smlc_1, smlc_2)

# label names are not given by a user
smlc_3 <- c(smlc1,smc1,smc1)
sml_IP4 <- ScoreMatrixList(list(IP1 = IP, IP2 = IP, IP, IP))
sml_control4 <- ScoreMatrixList(list(control, control, control, control))
smlc_4 <- ScoreMatrixListControl(sml_IP4, sml_control4 )
checkEquals(smlc_3, smlc_4)

# lack of the label name of one of the labels
smlc_5 <- c(smlc1, smc1, s2 = smc1)
sml_IP6 <- ScoreMatrixList(list(IP1 = IP, IP2 = IP, IP, s2 = IP))
sml_control6 <- ScoreMatrixList(list(control, control, control, control))
smlc_6 <- ScoreMatrixListControl(sml_IP6, sml_control6)
checkEquals(smlc_5, smlc_6)

# combine two ScoreMatrixListControl objects
smlc_7 <- c(smlc1, smlc1)
sml_IP8 <- ScoreMatrixList(list(IP1 = IP, IP2 = IP, IP1 = IP, IP2 = IP))
sml_control8 <- ScoreMatrixList(list(control, control, control, control))
smlc_8 <- ScoreMatrixListControl(sml_IP8, sml_control8)
checkEquals(smlc_7, smlc_8)

# combine a ScoreMatrixControl into a ScoreMatrixListControl object
smlc_9 <- c(smc1, smlc1)
sml_IP10 <- ScoreMatrixList(list(IP, IP1 = IP, IP2 = IP))
sml_control10 <- ScoreMatrixList(list(control, control, control))
smlc_10 <- ScoreMatrixListControl(sml_IP10, sml_control10)
checkEquals(smlc_9, smlc_10)

# combine ScoreMatrixControl objects into a ScoreMatrixListControl object
smlc_10 <- c(smc1, smlc1, smc1)
sml_IP11 <- ScoreMatrixList(list(IP, IP1 = IP, IP2 = IP, IP))
sml_control11 <- ScoreMatrixList(list(control, control, control, control))
smlc_11 <- ScoreMatrixListControl(sml_IP11, sml_control11)
checkEquals(smlc_10, smlc_11)

# combine ScoreMatrixControl objects into a ScoreMatrixListControl object
smlc_12 <- c(s1 = smc1, smlc1, s2 = smc1)
sml_IP13 <- ScoreMatrixList(list(s1 = IP, IP1 = IP, IP2 = IP, s2 = IP))
sml_control13 <- ScoreMatrixList(list(control, control, control, control))
smlc_13 <- ScoreMatrixListControl(sml_IP13, sml_control13)
checkEquals(smlc_12, smlc_13)

# ---------------------------------------------------------------------------- #
# enrichmentMatrix()
# ---------------------------------------------------------------------------- #

# ScoreMatrixControl
esm_1 <- enrichmentMatrix(smc1)
esm_2 <-log2((IP + 1) / (control + 1))
checkEquals(esm_1, esm_2)

# ScoreMatrixListControl
esml_1 <- enrichmentMatrix(smlc1)
es1 <- log2((IP + 1) / (control + 1))
es2 <- log2((IP + 1) / (control + 1))
esml_2 <- ScoreMatrixList(list(IP1 = es1, IP2 = es1))
checkEquals(esml_1, esml_2)

