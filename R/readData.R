# ------------------------------------------------------------------------------------ #

#' @param gl a \code{GRangesList} object, containing ranges for which represent regions enriched for transcription factor binding
#' @param width \code(integer) is the requested width of each enriched region. If 0 the ranges are not resized, if a positive integer, the width of all ranges is set to that number. Ranges are resized relative to the center of original ranges.
#' @param use.names a boolean which tells the function wheether to return the resulting ranges with a numeric vector which designates each class (the default), or to construct the names of each class using the names from the GRangesList
#' @param collapse.char a character which will be used to separate the class names if use.names=TRUE. The default is ':'

#' @usage readBroadPeak(path)
#' @return nothing

#' @docType methods
#' @rdname readBroadPeak
#' @export
setGeneric("readBroadPeak", function(path) standardGeneric("readBroadPeak") )

#' @aliases readBroadPeak
#' @rdname readBroadPeak
setMethod("readBroadPeak", signature("character"),
          function(path){
          
            # checks whether the file contains a track header
            header = scan(path, nmax=1, what='character', sep='\n')
            skip=0
            if(grepl('^track', header))
              skip=TRUE
            
            df = .readTableFast(path, skip=skip, header=FALSE, sep='\t')
            colnames(df) = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pvalue','qvalue')
            df$strand[df$strand == '.'] = '*'
            g = makeGRangesFromDataFrame(df=df, 
                                         keep.extra.columns=TRUE,
                                         seqnames.field = 'chrom',
                                         start.field = 'chromStart',
                                         end.field = 'chromEnd', 
                                         starts.in.df.are.0based=TRUE)
            return(g)
          }
          
          