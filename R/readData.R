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
#' @param remove.unsual if TRUE(default) remove the chromosomes with unsual 
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
#'  my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
#'  refseq = readGeneric(my.file,chr=1,start=2,end=3,strand=NULL,
#'                       meta.cols=list(score=5,name=4),
#'                      keep.all.metadata=FALSE, zero.based=TRUE)
#'  head(refseq)
#'
#' @export
#' @docType methods
#' @rdname readGeneric
readGeneric<-function(file, chr=1,start=2,end=3,strand=NULL,meta.col=NULL, 
                      keep.all.metadata=FALSE, zero.based=FALSE, 
                      remove.unsual=FALSE, header=FALSE, 
                      skip=0, sep="\t"){

  # removes the track header
  track=scan(file, n=1, what='character', sep='\n', quiet=TRUE)                    
  if(grepl('^>',track))
    skip = max(1, skip)
  
  # if header = FALSE, checks whether a header line exists, and skips it
  if(header==FALSE){
    line = read.table(file, nrows=1, stringsAsFactors=FALSE)
    if(all(sapply(line, class) == 'character'))
      skip = max(1, skip)
  }
  
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
  if(remove.unsual)
    df = df[grep("_", as.character(df$chr),invert=T),]
  
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
    mcols(g) = my.mcols[, !colnames(my.mcols) %in% black.names]
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
#' @examples
#' 
#' broad.peak.file = system.file('extdata',"ex.broadPeak", package='genomation')
#'                            
#' broad.peak = readBroadPeak(broad.peak.file)
#' head(broad.peak)
#'
#' @docType methods
#' @rdname readBroadPeak
#' @export
readBroadPeak<-function(file){
          
            g = readGeneric(file,
                            strand=6,
                            meta.cols=list(name=4,
                                          score=5,
                                          signalValue=7,
                                          pvalue=8,
                                          qvalue=9),
                            header=FALSE )
            return(g)
          }
        
# ---------------------------------------------------------------------------- #
#' A function to read the Encode formatted narrowPeak file into a GRanges object
#' @param file a abosulte or relative path to a bed file formatted by the Encode 
#'             narrowPeak standard
#' @usage readNarrowPeak(file)
#' @return a GRanges object
#'
#' @examples
#' 
#' narrow.peak.file = system.file('extdata',"ex.narrowPeak", package='genomation')
#'                  
#' narrow.peak = readBroadPeak(narrow.peak.file)
#' head(narrow.peak)
#'
#' @docType methods
#' @rdname readNarrowPeak
#' @export
readNarrowPeak<-function(file){
                      
            g = readGeneric(file,
                            strand=6,
                            meta.cols=list(name=4,
                                          score=5,
                                          signalValue=7,
                                          pvalue=8,
                                          qvalue=9,
                                          peak=10),
                            header=FALSE )
            return(g)
          }
         


# ---------------------------------------------------------------------------- #
#' A function to read-in genomic features and their upstream and downstream adjecent regions such as CpG islands and their shores
#'
#' @param location for the bed file of the feature 
#' @param flank number of basepairs for the flanking regions
#' @param clean If set to TRUE, flanks overlapping with other main features will be trimmed
#' @param remove.unsual  remove chromsomes with unsual names random, Un and antyhing with "_" character
#' @param feature.flank.name the names for feature and flank ranges, it should be a character vector of length 2. example: c("CpGi","shores")
#' @usage  readFeatureFlank(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL)
#' @return a GenomicRangesList contatining one GRanges object for flanks and one for GRanges object for the main feature.
#'   NOTE:This can not return a GRangesList at the moment because flanking regions do not have to have the same column name as the feature.
#'   GRangesList elements should resemble eachother in the column content. We can not satisfy that criteria for the flanks
#'
#' @examples
#'  cgi.path = system.file('extdata/chr21.CpGi.hg19.bed', package='genomation')
#'  cgi.shores = readFeatureFlank(cgi.path)
#'
#' @export
#' @docType methods
#' @rdname readFeatureFlank-methods
setGeneric("readFeatureFlank", 
           function(location,remove.unsual=T,
                    flank=2000,
                    clean=T,
                    feature.flank.name=NULL) 
             standardGeneric("readFeatureFlank") )

#' @aliases readFeatureFlank,character-method
#' @rdname readFeatureFlank-methods
setMethod("readFeatureFlank", 
          signature(location = "character"),
          function(location,remove.unsual,flank ,clean,feature.flank.name){
            
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
#' @param location location of the bed file with 12 or more columns
#' @param remove.unsual remove the chromomesomes with unsual names, mainly random chromsomes etc
#' @param up.flank  up-stream from TSS to detect promoter boundaries
#' @param down.flank down-stream from TSS to detect promoter boundaries
#' @param unique.prom     get only the unique promoters, promoter boundaries will not have a gene name if you set this option to be TRUE
#' @usage readTranscriptFeatures(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE)
#' @return a \code{\link{GRangesList}} containing locations of exon/intron/promoter/TSS
#' @note  one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#' 
#' @examples
#'   my.bed12.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
#'   my.bed12.file
#'   feats = readTranscriptFeatures(my.bed12.file) 
#'   names(feats)
#'   sapply(feats, head)
#' @export
#' @docType methods
#' @rdname readTranscriptFeatures-methods
setGeneric("readTranscriptFeatures", 
           function(location,remove.unsual=TRUE,
                    up.flank=1000,
                    down.flank=1000,
                    unique.prom=TRUE) 
             standardGeneric("readTranscriptFeatures"))

#' @aliases readTranscriptFeatures,character-method
#' @rdname readTranscriptFeatures-methods
setMethod("readTranscriptFeatures", 
          signature(location = "character"),
          function(location,remove.unsual,up.flank ,down.flank ,unique.prom){
            
            # find out if there is a header, skip 1st line if there is a header
            f.line=readLines(con = location, n = 1)
            skip=0
            if(grepl("^track",f.line))
              skip=1
            
            # readBed6
            message('Reading the table...\r')
            bed=readTableFast(location,header=F,skip=skip)                    
            if(remove.unsual)
              bed=bed[grep("_", as.character(bed[,1]),invert=T),]
            
            # introns
            message('Calculating intron coordinates...\r')
            introns	= convertBed2Introns(bed)
            # exons
            message('Calculating exon coordinates...\r')
            exons	= convertBed2Exons(bed)
            
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
#' @param gff.file path to a gff formatted file
#' @param split.group boolean, whether to split the 9th column of the file
#' @param split.char character that is used as a separator of the 9th column. ';' by default
#' @param filter a character designating which elements to retain from the gff file (e.g. exon, CDS, ...)
#' @param zero.based \Rcode{boolean} whether the coordinates are 0 or 1 based. 0 is the default
#' @return returns a \Rcode{GenomicRanges} object
#' 
#' @examples
#' gff.file = system.file('extdata/chr21.refseq.hg19.gtf', package='genomation')
#' gff = gffToGRanges(gff.file, split.group=TRUE)
#' 
#' 
#' @docType methods
#' @export
gffToGRanges = function(gff.file, split.group=FALSE, split.char=';',filter=NULL, 
                        zero.based=FALSE){
  
  gff = readGeneric(gff.file, 
                    chr=1,
                    start=4,
                    end=5, 
                    strand=7,
                    meta.cols=list(source=2,
                                  feature=3,
                                  score=6,
                                  frame=8,
                                  group=9), 
                    zero.based=zero.based)
  
  if(split.group){
    message('splitting the group.column...')
    group = strsplit(gff$group, '\\s+')
    gnames = group[[1]][seq(1,length(group[[1]]),2)]
    gids = seq(2,length(group[[1]]),2)
    group = data.frame(do.call(rbind, 
                              lapply(group, 
                                     function(x)gsub(split.char,'',x[gids]))), 
                       stringsAsFactors=FALSE)
    colnames(group) = gnames
    gff$group = NULL
    values(gff) = cbind(values(gff), group)
  }
  
  if(!is.null(filter)){		
    if(filter %in% gff$feature){
      message("Filtering", filter, "features...\n")
      gff = gff[gff$feature == filter,]
    }else{
      stop("The given feature is not present in the gff file")
    }
  }
  
  return(gff)
}

