test_annotateWithFeature = function(){
    target = GRanges(
        rep('chr1', 5),
        IRanges(c(1, 1, 6, 11, 16), c(1, 5, 10, 15, 20))
    )
    promoter = GRanges('chr1', IRanges(1, 5))
    intron = GRanges('chr1', IRanges(5, 10))
    exon = GRanges('chr1', IRanges(10, 15))

    res = annotateWithFeature(target, promoter)
    RUnit::checkIdentical(class(res)[1], "AnnotationByFeature")
    RUnit::checkIdentical(res@members, matrix(c(1, 1, 0, 0, 0), ncol = 1))
    RUnit::checkIdentical(res@annotation, c("promoter" = 40, "other" = 60))

    res = annotateWithFeature(target, intron)
    RUnit::checkIdentical(class(res)[1], "AnnotationByFeature")
    RUnit::checkIdentical(res@members, matrix(c(0, 1, 1, 0, 0), ncol = 1))
    RUnit::checkIdentical(res@annotation, c("intron" = 40, "other" = 60))

    res = annotateWithFeature(target, exon)
    RUnit::checkIdentical(class(res)[1], "AnnotationByFeature")
    RUnit::checkIdentical(res@members, matrix(c(0, 0, 1, 1, 0), ncol = 1))
    RUnit::checkIdentical(res@annotation, c("exon" = 40, "other" = 60))

    # test annotateWithFeatures
    feature = GRangesList("promotor" = promoter, "intron" =  intron, "exon" =  exon)
    mat = matrix(c(
        1, 1, 0, 0, 0,
        0, 1, 1, 0, 0,
        0, 0, 1, 1, 0
        ), ncol = 3
    )
    colnames(mat) = c("promotor", "intron", "exon")

    res = annotateWithFeatures(target, feature)
    RUnit::checkIdentical(class(res)[1], "AnnotationByFeature")
    RUnit::checkIdentical(res@members, mat)
    RUnit::checkIdentical(res@annotation, c("promotor" = 40, "intron" = 40, "exon" = 40))
}
