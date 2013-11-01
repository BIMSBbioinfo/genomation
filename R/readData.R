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
#' 
#' The function reads a tabular  text file that contains location and other information
#' on genomic features and returns a \code{\link{GRanges}} object. 
#' The minimal information that the file has to have is chromosome, 
#' start and end columns. Strand information is not compulsory.
#` 
#'  
#' @param file location of the file, a character string such as: "/home/user/my.bed"
#'
#' @param chr  number of the column that has chromsomes information in the table (Def:1)
#' @param start number of the column that has start coordinates in the table (Def:2)
#' @param end  number of the column that has end coordinates in the table (Def:3)
#' @param strand number of the column that has strand information, only -/+
#'               is accepted (Default:NULL)
#'               
#' @param meta.col  named \code{list} that maps column numbers to 
#'                   meta data columns. 
#'                   e.g. list(name=5, score=10), which means 5th column will be
#'                   named "name", and 10th column will be named "score" and their
#'                   contents will be a part of the returned GRanges object. 
#'                   If header = TRUE, meta.col parameter will over-write the 
#'                   column names given by the header line of the data frame.
#'                   
#' @param keep.all.metadata \code{logical} determining if the extra columns (
#'        the ones that are not designated by chr,start,end,strand and meta.col
#'        arguments )
#'        should be kept or not. (Default:FALSE)
#'        
#'        
#' @param zero.based a boolean which tells whether the ranges in 
#'        the bed file are 0 or 1 base encoded. (Default: FALSE)
#'        
#' @param remove.unsual if TRUE(default) remove the chromosomes with unsual 
#'        names, such as chrX_random (Default:FALSE)
#'        
#' @param header whether the original file contains a header line
#'  which designates the column names. If \code{TRUE} header will be used to 
#'  construct column names. These names can be over written by meta.col argument.
#' @param skip number of lines to skip. If there is a header line(s) you do not
#'        wish to include you can use skip argument to skip that line.
#' @param sep a single character which designates the separator in the file. 
#'        The default value is tab. 
#'        
#' @return \code{\link{GRanges}} object
#'
#' @examples
#'  my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
#'  readGeneric(my.file,chr=1,start=2,end=3,strand=NULL,
#'                       meta.col=list(score=5,name=4),
#'                      keep.all.metadata=FALSE, zero.based=TRUE)
#'
#' @export
#' @docType methods
#' @rdname readGeneric
readGeneric<-function(file, chr=1,start=2,end=3,strand=NULL,meta.col=NULL, 
                      keep.all.metadata=FALSE, zero.based=FALSE, remove.unsual=FALSE,
                      header=FALSE, skip=0,sep="\t"){
              
  # reads the bed files
  df=readTableFast(file, header=header, skip=skip, sep=sep)                    
  
  # removes nonstandard chromosome names
  if(remove.unsual)
    df = df[grep("_", as.character(df[,1]),invert=T),]
  
  # make a list of new column names, and their column numbers
  col.names1=list(chr=chr,start=start,end=end,strand=strand)
  col.names=c(col.names1,meta.col) # put the meta colums if any
  
  # check if col number exceeds dimensions of the original df.
  if( max(unlist(col.names)) > ncol(df) ) 
    stop("Number of columns is lower than designated number of columns by ",
         "meta.col,chr,start,end or strand arguments\n")
  
  # change the col names to the ones given by meta.col and chr,str,end,strand
  colnames(df)[unlist(col.names)] = names(unlist(col.names))
  
  
  g = makeGRangesFromDataFrame(
                          df, 
                          keep.extra.columns=FALSE, 
                          starts.in.df.are.0based=zero.based,
                          ignore.strand=is.null(strand))
  if(keep.all.metadata){
    mcols(g)=df[,-unlist(col.names1),drop=FALSE]
  }else if(!is.null(meta.col)){
    mcols(g)=df[,unlist(meta.col),drop=FALSE]
  }
    
  return(g)
}









# ---------------------------------------------------------------------------- #
#' A function to read the Encode formatted broad peak file into a GRanges object
#'
#' @param file a abosulte or relative path to a bed file formatted by the Encode broadPeak standard
#' @usage readBroadPeak(file)
#' @return a GRanges object
#'
#' @docType methods
#' @rdname readBroadPeak
#' @export
readBroadPeak<-function(file){
          
            # checks whether the file contains a track header
            #cnames = c('chrom','chromStart','chromEnd','name','score',
            #             'strand','signalValue','pvalue','qvalue')
            g = readGeneric(file,strand=6,meta.col=list(name=4,score=5,
                                                        signalValue=7,
                                                        pvalue=8,qvalue=9),
                            header=FALSE )
            return(g)
          }
        
# ---------------------------------------------------------------------------- #
#' A function to read the Encode formatted narrowPeak file into a GRanges object
#' @param file a abosulte or relative path to a bed file formatted by the Encode narrowPeak standard
#' @usage readNarrowPeak(file)
#' @return a GRanges object
#'        
#' @docType methods
#' @rdname readNarrowPeak
#' @export
readNarrowPeak<-function(file){
                      
            # checks whether the file contains a track header
            cnames = c('chrom','chromStart','chromEnd','name','score',
                         'strand','signalValue','pvalue','qvalue', 'peak')
            col.names = as.list(1:length(cnames))
            names(col.names) = cnames
            g = readGeneric(file,strand=6,meta.col=list(name=4,score=5,
                                                        signalValue=7,
                                                        pvalue=8,qvalue=9,
                                                        peak=10),
                            header=FALSE )
            return(g)
          }
         
          