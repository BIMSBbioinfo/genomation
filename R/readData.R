# ---------------------------------------------------------------------------- #
# fast reading of big tables
readTableFast<-function(filename,header=T,skip=0,sep=""){
  
  tab5rows <- read.table(filename, header = header,skip=skip,
                         sep=sep, nrows = 100, stringsAsFactors=F)
  classes  <- sapply(tab5rows, class)
  df = read.table(filename, 
                  header = header,
                  skip=skip,
                  sep=sep, 
                  colClasses = classes,
                  stringsAsFactors=FALSE)
  return(df)
}

# ---------------------------------------------------------------------------- #
#' Read a tabular file and convert it to GRanges. 
#' The function can take a bed file with or without a predesignated header, or 
#' any other tabular file.
#' If the file is in bed format, the user must specifiy which columns contain
#' chromosome, start, and end information, using the col.names argument. 
#' Strand information is not compulsory.
#` 
#'  
#' @param file location of the file, a character string such as: "/home/user/my.bed"
#' @param header whether the original file contains a header line which designates the column names. 
#' @param col.names a named \Rcode{list} that maps individual columns to chr, start, end, strand information. e.g. list(chr=4, start=6, end=7, strand=9, score=10). If header = TRUE, colnames parameter does not have any effect. If the keep.medata is set to FALSE, only chromosome, start, end and strand colums will be used to construct the GRanges object, otherwise, all columns will be read in.
#' 
#' @param starts.in.df.are.0based a boolean which tells whether the ranges in the bed file are 0 or 1 base encoded. A 0 based encoding is persumed
#' @param remove.unsual if TRUE(default) remove the chromomesomes with unsual names, such as chrX_random
#' @param sep a single character which designates the separator in the file. The default value is tab. 
#'
#' @usage readGeneric(location,remove.unsual=T)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname readGeneric-methods
setGeneric("readGeneric", 
            function(file, 
                     header=FALSE, 
                     col.names=NULL,
                     keep.metadata=TRUE,
                     starts.in.df.are.0based=TRUE,
                     remove.unsual=TRUE,
                     sep='\t',
                     ...) 
              standardGeneric("readGeneric"))

#' @aliases readGeneric,character-method
#' @rdname readGeneric-methods
setMethod("readGeneric", 
          signature(file = "character"),
          function(file, header, col.names, keep.metadata, 
                   starts.in.df.are.0.based, remove.unsual, sep, ...){
            
            # find out if there is a header, skip 1st line if there is a header
            f.line=readLines(con = file, n = 1)
            skip=0
            if(grepl("^track",f.line))
              skip=1
            
            # reads the bed files
            bed=readTableFast(file, header=header, skip=skip, sep=sep)                    
            
            # removes nonstandard chromosome names
            if(remove.unsual)
              bed = bed[grep("_", as.character(bed[,1]),invert=T),]
            
            # checks whether the bed file contains a designated header, 
            # if it does not, then it uses the col.names variable to construct the header
            # if col.names is empty it will name the first three columns chr start end
            if(!header){
              if(is.null(col.names)){
                colnames(bed)[1:3] = c('chr','start','end')
              }else{
                if(!is.list(col.names))
                  stop('col.names argument is not a named list')
                  
                  # checks whether all necessary information is present
                  if(!all(c('chr','start','end') %in% names(col.names)))
                    stop(paste(setdiff(c('chr','start','end'), names(col.names))
                               ,'information is missing' ))
                  
                  colnames(bed)[unlist(col.names)] = names(col.names)
                }
              }
            
            
            # converts . as strand character into *
            sind = grepl('strand', colnames(bed))
            if(sind)
              bed[, sind] = gsub('\\.','*', bed[,sind])
            
            g = makeGRangesFromDataFrame(
                                    bed, 
                                    keep.extra.columns=keep.metadata, 
                                    starts.in.df.are.0based=starts.in.df.are.0based,
                                    ignore.strand=!sind)
            return(g)
          }
)


# ---------------------------------------------------------------------------- #
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
            cnames = c('chrom','chromStart','chromEnd','name','score',
                         'strand','signalValue','pvalue','qvalue')
            col.names = as.list(1:length(cnames))
            names(col.names) = cnames
            
            g = readGeneric(file, header=FALSE, col.names=col.names)
            return(g)
          }
)          
# ---------------------------------------------------------------------------- #
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
            cnames = c('chrom','chromStart','chromEnd','name','score',
                         'strand','signalValue','pvalue','qvalue', 'peak')
            col.names = as.list(1:length(cnames))
            names(col.names) = cnames
            g = readGeneric(file, header=FALSE, colnames=colnames)
            return(g)
          }
)          
          