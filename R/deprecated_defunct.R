#' @export
plotGeneAnnotation  <- function(l, 
                                cluster=FALSE, 
                                col=c('white','cornflowerblue'))
{
  .Deprecated("heatTargetAnnotation ")
  ## use new function, or remainder of myOldFunc
  heatTargetAnnotation(l,cluster=cluster, col=col,
                         precedence=FALSE,plot=TRUE)
}