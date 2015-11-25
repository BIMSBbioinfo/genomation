# ---------------------------------------------------------------------------- #

read.zip <- function(file, ...) {
  zipFileInfo <- unzip(file, list=TRUE)
  if(nrow(zipFileInfo) > 1)
    stop("More than one data file inside zip")
  else
    read.table(unz(file, as.character(zipFileInfo$Name)), ...)
}

compressedAndUrl2temp <- function(filename){ 
  if(grepl('^(http://|https://|ftp://|ftps://).*(.gz|.bz2|.xz|.zip)$',filename)){
    temp <- tempfile()
    download.file(filename, destfile=temp, method="auto", quiet = TRUE)
    ext <- tail(unlist(strsplit(x=filename, split="\\.")), n=1)
    filename <- paste0(temp, ".", ext)
    file.rename(temp, filename)
  }
  filename
}
# detects UCSC header (and first track)
detectUCSCheader <- function(filename){
  skip=0
  if(grepl("^.*(.zip)[[:space:]]*$", filename)){
    zipFileInfo <- unzip(filename, list=TRUE)
    if(length(zipFileInfo$Name)>1) stop("More than one data file inside zip")
    con <- unz(filename, filename=zipFileInfo$Name); open(con)
  }else{
    con <- file(filename, open = "r")
  }
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(grepl("^ *browser", oneLine) | 
         grepl("^ *track", oneLine) |
         grepl("^ *#", oneLine)){
      skip = skip + 1
    }else{
      break
    }
  }  
  close(con)
  return(skip)
}


# fast reading of big tables
readTableFast<-function(filename,header=TRUE,skip=0,sep="\t"){
  
  filename <- compressedAndUrl2temp(filename)
  if(skip==FALSE){
    skip = 0
  }else if(skip=="auto"){
    skip = detectUCSCheader(filename)
  }
 
  if(grepl("^.*(.zip)[[:space:]]*$", filename)){
	tab30rows <- read.zip(filename,
			      sep=sep, 
			      skip=skip,
			      nrows=30,
			      header=header,
			      stringsAsFactors=FALSE)
  }else{
	tab30rows <- read.table(file=filename, 
				sep=sep, 
				skip=skip,
				nrows=30,
				header=header,
				stringsAsFactors=FALSE)
  }
  classes  <- sapply(tab30rows, class)
  cl=""
  for(cla in classes){
    if(cla=="character"){
      cl=paste0(cl, "c")
    }else if(cla=="numeric" | cla=="double"){
      cl=paste0(cl, "d")
    }else if(cla=="integer"){
      cl=paste0(cl, "i")
    }else if(cla=="logical")
      cl=paste0(cl, "l")
  }
  
  df <- read_delim(file=filename, 
                   delim=sep, 
                   skip=skip,
                   col_names=header,
                   col_types=cl,
                   locale=locale(grouping_mark = "_")
                   )
  #changing default variables names from read_delim (X[0-9]+) to data.frame (V[0-9]+)
  colnames(df) <- gsub("^X(\\d+)$", "V\\1", colnames(df)) 
  return(as.data.frame(df))
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
#'             or the input itself as a string (containing at least one \\n).
#'             
#' @param chr  number of the column that has chromsomes information in the table (Def:1)
#' @param start number of the column that has start coordinates in the table (Def:2)
#' @param end  number of the column that has end coordinates in the table (Def:3)
#' @param strand number of the column that has strand information, only -/+
#'               is accepted (Default:NULL)
#'               
#' @param meta.cols  named \code{list} that maps column numbers to 
#'                   meta data columns. 
#'                   e.g. list(name=5, score=10), which means 5th column will be
#'                   named "name", and 10th column will be named "score" and their
#'                   contents will be a part of the returned GRanges object. 
#'                   If header = TRUE, meta.cols parameter will over-write the 
#'                   column names given by the header line of the data frame.
#'                   
#' @param keep.all.metadata \code{logical} determining if the extra columns (
#'        the ones that are not designated by chr,start,end,strand and meta.cols
#'        arguments )
#'        should be kept or not. (Default:FALSE)
#'        
#'        
#' @param zero.based a boolean which tells whether the ranges in 
#'        the bed file are 0 or 1 base encoded. (Default: FALSE)
#'        
#' @param remove.unusual if TRUE(default) remove the chromosomes with unsual 
#'        names, such as chrX_random (Default:FALSE)
#'        
#' @param header whether the original file contains a header line
#'  which designates the column names. If \code{TRUE} header will be used to 
#'  construct column names. These names can be over written by meta.cols argument.
#' @param skip number of lines to skip. If there is a header line(s) you do not
#'        wish to include you can use skip argument to skip that line.
#' @param sep a single character which designates the separator in the file. 
#'        The default value is tab. 
#'        
#' @return \code{\link{GRanges}} object
#'
#' @examples
#' my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
#' refseq = readGeneric(my.file,chr=1,start=2,end=3,strand=NULL,
#'                       meta.cols=list(score=5,name=4),
#'                       keep.all.metadata=FALSE, zero.based=TRUE)
#' head(refseq)
#' @export
#' @docType methods
#' @rdname readGeneric
readGeneric<-function(file, chr=1,start=2,end=3,strand=NULL,meta.cols=NULL, 
                      keep.all.metadata=FALSE, zero.based=FALSE, 
                      remove.unusual=FALSE, header=FALSE, 
                      skip=0, sep="\t"){
  
  # reads the bed files
  df=readTableFast(file, header=header, skip=skip, sep=sep)                    
  
  # make a list of new column names, and their column numbers
  col.names1=list(chr=chr,start=start,end=end,strand=strand)
  col.names=c(col.names1,meta.cols) # put the meta colums if any
  
  # check if col number exceeds dimensions of the original df.
  if( max(unlist(col.names)) > ncol(df) ) 
    stop("Number of columns is lower than designated number of columns by ",
         "meta.cols,chr,start,end or strand arguments\n")
  
  # change the col names to the ones given by meta.cols and chr,str,end,strand
  colnames(df)[unlist(col.names)] = names(unlist(col.names))
  
  # converts the . strand character to *
  sind = grepl('strand',colnames(df))
  if(any(sind) & !is.null(strand))
    df[, sind] = sub('\\.', '*', df[,sind])
  
  # removes nonstandard chromosome names
  if(remove.unusual)
    df = df[grep("_", as.character(df$chr),invert=TRUE),]
  
  g = makeGRangesFromDataFrame(
    df, 
    keep.extra.columns=FALSE, 
    starts.in.df.are.0based=zero.based,
    ignore.strand=is.null(strand))
  
  # this names can not be column names in meta-data
  black.names=c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                "isCircular", "start", "end", "width", "element")
  
  
  if(keep.all.metadata){
    my.mcols = df[,-unlist(col.names1),drop=FALSE]
    mcols(g) = my.mcols[, !colnames(my.mcols) %in% black.names]
  }else if(!is.null(meta.cols)){
    my.mcols=df[,unlist(meta.cols),drop=FALSE]
    values(g) = my.mcols[, !colnames(my.mcols) %in% black.names, drop=FALSE]
  }
  
  return(g)
}

#' Read a BED file and convert it to GRanges. 
#' 
#' The function reads a BED file that contains location and other information
#' on genomic features and returns a \code{\link{GRanges}} object. 
#' The minimal information that the BED file has to have is chromosome, 
#' start and end columns. it can handle all BED formats up to 12 columns. 
#' 
#'  
#' @param file location of the file, a character string such as: "/home/user/my.bed"
#'             or the input itself as a string (containing at least one \\n).
#'             The file can end in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip}
#'             and/or start with \code{http://} or \code{ftp://}. If the file is not 
#'             compressed it can also start with \code{https://} or \code{ftps://}.
#'
#' @param track.line the number of track lines to skip, "auto" to detect them 
#'                   automatically or FALSE(default) if the bed file doesn't 
#'                   have track lines
#'                   
#' @param remove.unusual if TRUE remove the chromosomes with unsual 
#'        names, such as chrX_random (Default:FALSE)        
#' @param zero.based a boolean which tells whether the ranges in 
#'        the bed file are 0 or 1 base encoded. (Default: TRUE)
#'        
#' @return \code{\link{GRanges}} object
#'
#' @examples
#' my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
#' refseq = readBed(my.file,track.line=FALSE,remove.unusual=FALSE)
#' head(refseq)
#' 
#' @export
#' @docType methods
#' @rdname readBed
readBed<-function(file,track.line=FALSE,remove.unusual=FALSE,zero.based=TRUE)
{
  
  meta.cols=list(score=5,name=4,thickStart=7,  
                 thickEnd=8, 
                 itemRgb=9,  
                 blockCount=10, 
                 blockSizes=11,  
                 blockStarts=12 )
  
  file <- compressedAndUrl2temp(file)
  if(is.numeric(track.line)){
    skip = track.line
  }else if(track.line=="auto"){
    skip = detectUCSCheader(file)
  }else{
    skip = 0
  }
  
  df=read_delim(file,skip=skip,n_max=2,col_names=FALSE, delim="\t")
  numcol=ncol(df)
  
  if(numcol==3){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL,   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")
  }else if(numcol==4){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = meta.cols[1],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")    
  }else if(numcol==5){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = meta.cols[1:2],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")    
  }else if(numcol == 6){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = 6,
                   meta.cols = meta.cols[1:2],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t") 
  }else if(numcol > 6){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = 6,
                   meta.cols = meta.cols[c(1:2,(3:8)[1:(numcol-6)])],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t") 
  }
  df
  
} 

# ---------------------------------------------------------------------------- #
#' A function to read the Encode formatted broad peak file into a GRanges object
#'
#' @param file an absolute or relative path to a bed file formatted by the Encode broadPeak standard.
#'             The file can end in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip}
#'             and/or start with \code{http://} or \code{ftp://}. If the file is not compressed
#'             it can also start with \code{https://} or \code{ftps://}.
#' @param track.line the number of track lines to skip, "auto" to detect them automatically
#'                   or FALSE(default) if the bed file doesn't have track lines
#' @usage readBroadPeak(file, track.line=FALSE)
#' @return a GRanges object
#'
#' @examples
#' broad.peak.file = system.file('extdata',"ex.broadPeak", package='genomation')
#'                            
#' broad.peak = readBroadPeak(broad.peak.file)
#' head(broad.peak)
#' 
#' @docType methods
#' @rdname readBroadPeak
#' @export
readBroadPeak<-function(file, track.line=FALSE){
  
  g = readGeneric(file,
                  strand=6,
                  meta.cols=list(name=4,
                                 score=5,
                                 signalValue=7,
                                 pvalue=8,
                                 qvalue=9),
                  header=FALSE,
                  skip=track.line)
  return(g)
}

# ---------------------------------------------------------------------------- #
#' A function to read the Encode formatted narrowPeak file into a GRanges object
#' @param file an absolute or relative path to a bed file formatted by the Encode 
#'             narrowPeak standard. The file can end in \code{.gz}, \code{.bz2}, 
#'             \code{.xz}, or \code{.zip} and/or start with \code{http://} or 
#'             \code{ftp://}. If the file is not compressed
#'             it can also start with \code{https://} or \code{ftps://}.
#' @param track.line the number of track lines to skip, "auto" to detect them automatically
#'                   or FALSE(default) if the bed file doesn't have track lines
#' @usage readNarrowPeak(file, track.line=FALSE)
#' @return a GRanges object
#'
#' @examples
#' narrow.peak.file = system.file('extdata',"ex.narrowPeak", package='genomation')
#'                  
#' narrow.peak = readBroadPeak(narrow.peak.file)
#' head(narrow.peak)
#' @docType methods
#' @rdname readNarrowPeak
#' @export
readNarrowPeak<-function(file, track.line=FALSE){
  
  g = readGeneric(file,
                  strand=6,
                  meta.cols=list(name=4,
                                 score=5,
                                 signalValue=7,
                                 pvalue=8,
                                 qvalue=9,
                                 peak=10),
                  header=FALSE,
                  skip=track.line)
  return(g)
}



# ---------------------------------------------------------------------------- #
#' A function to read-in genomic features and their upstream and downstream adjecent regions
#' such as CpG islands and their shores
#'
#' @param location for the bed file of the feature.
#' @param flank number of basepairs for the flanking regions
#' @param clean If set to TRUE, flanks overlapping with other main features will be trimmed
#' @param remove.unusual  remove chromsomes with unsual names random, Un and antyhing with "_" character
#' @param feature.flank.name the names for feature and flank ranges, it should be a character 
#'                           vector of length 2. example: c("CpGi","shores")
#' @usage readFeatureFlank(location,remove.unusual=TRUE,flank=2000,
#'                         clean=TRUE,feature.flank.name=NULL)
#' @return a GenomicRangesList contatining one GRanges object for flanks and one for GRanges object
#'         for the main feature.
#'   NOTE: This can not return a GRangesList at the moment because flanking regions do not
#'   have to have the same column name as the feature. GRangesList elements should resemble 
#'   each other in the column content. We can not satisfy that criteria for the flanks
#'
#' @examples
#' cgi.path = system.file('extdata/chr21.CpGi.hg19.bed', package='genomation')
#' cgi.shores = readFeatureFlank(cgi.path)
#' cgi.shores
#' @export
#' @docType methods
#' @rdname readFeatureFlank-methods
setGeneric("readFeatureFlank", 
           function(location,remove.unusual=TRUE,
                    flank=2000,
                    clean=TRUE,
                    feature.flank.name=NULL) 
             standardGeneric("readFeatureFlank") )

#' @aliases readFeatureFlank,character-method
#' @rdname readFeatureFlank-methods
setMethod("readFeatureFlank", 
          signature(location = "character"),
          function(location,remove.unusual,flank ,clean,feature.flank.name){
            
            feat = readGeneric(location)
            flanks = getFlanks(feat,flank=flank,clean=clean)
            x = GenomicRangesList(features=feat,flanks=flanks)
            if(!is.null(feature.flank.name) & length(feature.flank.name)==2)
              names(x)=feature.flank.name
            
            return(x)
          })

# ---------------------------------------------------------------------------- #
#' Function for reading exon intron and promoter structure from a given bed file
#'
#' @param location location of the bed file with 12 or more columns. 
#'                 The file can end in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip}
#'                 and/or start with \code{http://} or \code{ftp://}. If the file is not compressed
#'                 it can also start with \code{https://} or \code{ftps://}.
#' @param remove.unusual remove the chromomesomes with unsual names, mainly random chromsomes etc
#' @param up.flank  up-stream from TSS to detect promoter boundaries
#' @param down.flank down-stream from TSS to detect promoter boundaries
#' @param unique.prom get only the unique promoters, promoter boundaries will not have 
#'                    a gene name if you set this option to be TRUE
#' @usage readTranscriptFeatures(location,remove.unusual=TRUE,
#'                               up.flank=1000,down.flank=1000,unique.prom=TRUE)
#' @return a \code{\link{GRangesList}} containing locations of exon/intron/promoter/TSS
#' @note  one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#' 
#' @examples
#' my.bed12.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
#' my.bed12.file
#' feats = readTranscriptFeatures(my.bed12.file) 
#' names(feats)
#' sapply(feats, head)
#' @export
#' @docType methods
#' @rdname readTranscriptFeatures-methods
setGeneric("readTranscriptFeatures", 
           function(location,remove.unusual=TRUE,
                    up.flank=1000,
                    down.flank=1000,
                    unique.prom=TRUE) 
             standardGeneric("readTranscriptFeatures"))

#' @aliases readTranscriptFeatures,character-method
#' @rdname readTranscriptFeatures-methods
setMethod("readTranscriptFeatures", 
          signature(location = "character"),
          function(location,remove.unusual,up.flank ,down.flank ,unique.prom){
            
            # readBed6
            message('Reading the table...\r')
            bed=readTableFast(location,header=FALSE,skip="auto")                    
            if(remove.unusual)
              bed=bed[grep("_", as.character(bed[,1]),invert=TRUE),]
            
            # introns
            message('Calculating intron coordinates...\r')
            introns    = convertBed2Introns(bed)
            # exons
            message('Calculating exon coordinates...\r')
            exons    = convertBed2Exons(bed)
            
            # get the locations of TSSes
            message('Calculating TSS coordinates...\r')
            tss=bed
            #  + strand
            tss[tss$V6=="+",3] = tss[tss$V6=="+",2]
            #  - strand
            tss[tss$V6=="-",2]=tss[tss$V6=="-",3]
            
            tssg = GRanges(seqnames=as.character(tss$V1),
                           ranges=IRanges(start=tss$V2, end=tss$V3),
                           strand=as.character(tss$V6),
                           score=rep(0,nrow(tss)),
                           name=tss$V4)
            
            message('Calculating promoter coordinates...\r')
            # get the locations of promoters
            # + strand
            bed[bed$V6=="+",3]=bed[bed$V6=="+",2]+down.flank
            bed[bed$V6=="+",2]=bed[bed$V6=="+",2]-up.flank
            
            #  - strand
            bed[bed$V6=="-",2]=bed[bed$V6=="-",3]-down.flank
            bed[bed$V6=="-",3]=bed[bed$V6=="-",3]+up.flank
            
            
            if(! unique.prom){
              prom.df = (bed[,c(1,2,3,4,6)])
              prom = GRanges(seqnames=as.character(prom.df$V1),
                             ranges = IRanges(start=prom.df$V2, end=prom.df$V3),
                             strand = as.character(prom.df$V6),
                             score=rep(0,nrow(prom.df)),
                             name=prom.df$V4)
            }else{
              prom.df = unique(bed[,c(1,2,3,6)])
              prom = GRanges(seqnames=as.character(prom.df$V1),
                             ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
                             strand=as.character(prom.df$V6),
                             score=rep(0,nrow(prom.df)),
                             name=rep(".",nrow(prom.df)) )
            }
            
            message('Outputting the final GRangesList...\r\n')
            GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
          })


# ---------------------------------------------------------------------------- #
#' Converts a gff formated data.frame into a GenomicRanges object. 
#' The GenomicRanges object needs to be properly formated for the function to work.
#' @param gff.file path to a gff formatted file. 
#'                 The file can end in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip}
#'                 and/or start with \code{http://} or \code{ftp://}. If the file is not compressed
#'                 it can also start with \code{https://} or \code{ftps://}.
#' @param filter a character designating which elements to retain from the gff file (e.g. exon, CDS, ...)
#' @param zero.based \code{boolean} whether the coordinates are 0 or 1 based. 0 is the default
#' @param ensembl \code{boolean} if TRUE, add the chr prefix to seqlevels. FALSE by default
#' @return returns a \code{GenomicRanges} object
#' 
#' @examples
#' gff.file = system.file('extdata/chr21.refseq.hg19.gtf', package='genomation')
#' gff = gffToGRanges(gff.file)
#' 
#' @docType methods
#' @export
gffToGRanges = function(gff.file, filter=NULL, zero.based=FALSE, ensembl=FALSE){
  
  gff = rtracklayer::import(gff.file)
  if(zero.based)
    gff$start = gff$start + 1
  
  if(ensembl)
    seqlevels(gff) = paste('chr',seqlevels(gff),sep='')
  
  if(!is.null(filter)){
    if(!any(gff$type == filter))
      stop(paste(filter, 'category does not exist in the gff file'))
    gff = gff[grepl(filter, gff$type)]
  }
  
  return(gff)
}
