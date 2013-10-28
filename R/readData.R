# ----------------------------------------------------------------------------------------------- #
# fast reading of big tables
.readTableFast<-function(filename, header=T, skip=0, sep=""){
  tab5rows <- read.table(filename, header = header,skip=skip,sep=sep, nrows = 100)
  classes  <- sapply(tab5rows, class)
  return( read.table(filename, header = header,skip=skip,sep=sep, colClasses = classes)  )
}


# ----------------------------------------------------------------------------------------------- #
#' read a bed file and convert it to GRanges
#'  
#' @param file  location of the file, a character string such as: "/home/user/my.bed"
#' @param remove.unsual if TRUE(default) remove the chromomesomes with unsual names, mainly random chromsomes etc
#'
#' @usage readBed(location,remove.unsual=T)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname readBed-methods
setGeneric("readBed", function(file,remove.unsual=T) standardGeneric("readBed"))

#' @aliases readBed,character-method
#' @rdname readBed-methods
setMethod("readBed", 
          signature(file = "character"),
          function(file, remove.unsual){
            
            # find out if there is a header, skip 1st line if there is a header
            f.line=readLines(con = file, n = 1)
            skip=0
            if(grepl("^track",f.line))
              skip=1
            
            # readBed6
            bed=.readTableFast(file,header=F,skip=skip)                    
            if(remove.unsual)
              bed = bed[grep("_", as.character(bed[,1]),invert=T),]
            
            convertBedDf(bed)
          })


# ----------------------------------------------------------------------------------------------- #
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
          
          