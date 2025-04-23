test_readTranscriptFeatures = function(){
    bed6 = system.file("unitTests/bed6.bed", package = "genomation")
    checkException(suppressMessages(readTranscriptFeatures(bed6)))

    bed12 = system.file("unitTests/bed12.bed", package = "genomation")
    features = suppressMessages(readTranscriptFeatures(bed12))

    checkTrue(inherits(features, "GenomicRangesList"))
    checkIdentical(c("exons", "introns", "promoters", "TSSes"), names(features))

    exons = as.data.frame(features$exons)
    exons_expected = data.frame(
        seqnames = factor("ch1"),
        start = c(12001L, 14001L, 16001L, 22001L, 24001L, 26001L),
        end = c(13000L, 15000L, 17000L, 23000L, 25000L, 27000L),
        width = 1000L,
        strand = factor(c("+", "+", "+", "-", "-", "-"), levels = c("+", "-", "*")),
        score = as.numeric(c(1:3, 3:1)),
        name = c(rep("f1", 3), rep("f2", 3))
    )
    checkIdentical(exons_expected, exons)

    introns = as.data.frame(features$introns)
    introns_expected = data.frame(
        seqnames = factor("ch1"),
        start = c(13001L, 15001L, 23001L, 25001L),
        end = c(14000L, 16000L, 24000L, 26000L),
        width = 1000L,
        strand = factor(c("+", "+", "-", "-"), levels = c("+", "-", "*")),
        score = c(1, 2, 2, 1),
        name = c("f1", "f1", "f2", "f2")
    )
    checkIdentical(introns_expected, introns)

    promoters = as.data.frame(features$promoters)
    promoters_expected = data.frame(
        seqnames = factor("ch1"),
        start = c(9000L, 29000L),
        end = c(11000L, 31000L),
        width = 2001L,
        strand = factor(c("+", "-"), levels = c("+", "-", "*")),
        score = 0,
        name = "."
    )
    checkIdentical(promoters_expected, promoters)

    TSSes = as.data.frame(features$TSSes)
    TSSes_expected = data.frame(
        seqnames = factor("ch1"),
        start = c(10000L, 30000L),
        end = c(10000L, 30000L),
        width = 1L,
        strand = factor(c("+", "-"), levels = c("+", "-", "*")),
        score = 0,
        name = c("f1", "f2")
    )
    checkIdentical(TSSes_expected, TSSes)
}
