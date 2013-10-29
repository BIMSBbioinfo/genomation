# ----------------------------------------------------------------------------------------------- #
# fast reading of big tables
readTableFast<-function(filename,header=T,skip=0,sep=""){
  
  tab5rows <- read.table(filename, header = header,skip=skip,sep=sep, nrows = 100, stringsAsFactors=F)
  classes  <- sapply(tab5rows, class)
  df = read.table(filename, 
                  header = header,
                  skip=skip,
                  sep=sep, 
                  colClasses = classes,
                  stringsAsFactors=FALSE)
  return(df)
}

# ----------------------------------------------------------------------------------------------- #
#' read a bed file and convert it to GRanges. 
#' The function can take a bed file with or without predesignated headers. If the file contains a track header line it will automatically skip it. 
#'  
#' @param file  location of the file, a character string such as: "/home/user/my.bed"
#' @param remove.unsual if TRUE(default) remove the chromomesomes with unsual names, mainly random chromsomes etc
#' @param colnames a character vector that designates the names of the columns. If empty string a standard ordering of columns is assumed (e.g. chromosome, start, end)
#' @param header a boolean whether the original file contains a header line which designates the column names.
#' @param colnames a character vector that gives designates the name of the columns of the bed file. Only applicable if the header is set to \Rcode{FALSE}. If header == FALSE and nchar(colnames) == 0, the function will set the names of the first three columns to chr, start, end
#' @param starts.in.df.are.0based a boolean which tells the whether the ranges in the bed file are 0 or 1 base encoded. A 0 based encoding is persumed
#'
#' @usage readBed(location,remove.unsual=T)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname readBed-methods
setGeneric("readBed", function(file, remove.unsual=TRUE, header=FALSE, colnames='', starts.in.df.are.0based=TRUE) standardGeneric("readBed"))

#' @aliases readBed,character-method
#' @rdname readBed-methods
setMethod("readBed", 
          signature(file = "character"),
          function(file, remove.unsual, header, colnames, starts.in.df.are.0.based){
            
            # find out if there is a header, skip 1st line if there is a header
            f.line=readLines(con = file, n = 1)
            skip=0
            if(grepl("^track",f.line))
              skip=1
            
            # reads the bed files
            bed=readTableFast(file, header=header, skip=skip)                    
            
            # removes nonstandard chromosome names
            if(remove.unsual)
              bed = bed[grep("_", as.character(bed[,1]),invert=T),]
            
            # checks whether the bed file contains a designated header, 
            # if it does not, then it uses the colnames variable to construct the header
            # if col.names is empty it will name the first three columns chr start end
            if(!header){
              if(nchar(col.names) == 0){
                colnames(bed)[1:3] = c('chr','start','end')
              }else{
                if(length(colnames(bed) != ncol(bed)))
                  stop('number of designated col.names does not equal the number of columns')
                colnames(bed)=col.names
              }
            }
            
            sind = grepl('strand', colnames(bed))
            if(sind)
              bed[, sind] = gsub('\\.','*', bed[,sind])
            
           g = makeGRangesFromDataFrame(bed, 
                                     keep.extra.columns=TRUE, 
                                     starts.in.df.are.0based=starts.in.df.are.0based,
                                     ignore.strand=!sind)
            return(g)
          }
)


# ----------------------------------------------------------------------------------------------- #
#` A function to read the Encode formatted broad peak file into a GRanges object
#' @param file a abosulte or relative path to a bed file formatted by the Encode broadPeak standard
#' @usage readBroadPeak(file)
#' @return a GRanges object

#' @docType methods
#' @rdname readBroadPeak
#' @export
setGeneric("readBroadPeak", function(file) standardGeneric("readBroadPeak") )

#' @aliases readBroadPeak
#' @rdname readBroadPeak
setMethod("readBroadPeak", signature("character"),
          function(file){
          
            # checks whether the file contains a track header
            colnames = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pvalue','qvalue')
            g = readBed(file, header=FALSE, colnames=colnames)
            return(g)
          }
)          
# ----------------------------------------------------------------------------------------------- #
#` A function to read the Encode formatted narrowPeak file into a GRanges object
#' @param file a abosulte or relative path to a bed file formatted by the Encode narrowPeak standard
#' @usage readNarrowPeak(file)
#' @return a GRanges object
          
#' @docType methods
#' @rdname readNarrowPeak
#' @export
setGeneric("readNarrowPeak", function(file) standardGeneric("readNarrowPeak") )
    
#' @aliases readNarrowPeak
#' @rdname readNarrowPeak
setMethod("readNarrowPeak", signature("character"),
          function(file){
                      
            # checks whether the file contains a track header
            colnames = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pvalue','qvalue', 'peak')
            g = readBed(file, header=FALSE, colnames=colnames)
            return(g)
          }
)          
          