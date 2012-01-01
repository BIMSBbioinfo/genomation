

#######################################
# S3 functions
#######################################

# returns a matrix
make.scoreMatrixBin<-function(target,windows,bin.num,bin.op)
{
            # chop windows to bins
            my.binner=function(strt,end,nbins){
              x=round(seq(from = strt, to = end,length.out=nbins ) )
              my.strt=x[1:(nbins-1)]
              my.end =x[2:nbins]-1
              
              return( t(cbind(my.strt, my.end) )  )
            }

            coord=matrix(
                        mapply(my.binner,start(windows),end(windows),bin.num+1,SIMPLIFY=TRUE )
                    ,ncol=2,byrow=T)
           subWins=GRanges(seqnames=rep(as.character(seqnames(windows)),bin.num),IRanges(start=coord[,1],end=coord[,2]))
            
            #check if there are common chromsomes
            win.list=as(subWins, "RangesList")
            chrs.t=names(target)
            chrs.w=names(win.list)
            chrs  =chrs.t[chrs.t %in% chrs.w]
            if(length(chrs)==0){
              stop("There are no common chromosomes/spaces to do overlap")
            }
            
            
            #get views           
            my.vList=Views( target[names(target) %in% chrs], win.list[names(win.list) %in% chrs] )            
            #RleViewsList( rleList = target[names(target) %in% chrs], rangesList = win.list[names(win.list) %in% chrs] )

            if(is(target,"modRleList") ){
              # get vectors from Views and make a matrix outof it

              #stores the number of indices with values in Rle vector list
              my.lList=Views( target[names(target) %in% chrs]>0, win.list[names(win.list) %in% chrs] )
              
              if(bin.op=="mean")
              {
                sum.bins=unlist(lapply(my.vList,viewSums) )# sum of each view
                len.bins=unlist(lapply(my.lList,viewSums ) ) # number of values in each bin, discarding bases with no value
                mat=matrix( ((sum.bins/len.bins)-target@add)/target@multiply, ncol=bin.num,byrow=TRUE)
              
                return( mat )
              }
              else if (bin.op=="max")
              {
                sum.bins=unlist(lapply(my.vList,viewMaxs) )# max of each view
              }
              else if (bin.op=="min")
              {
                sum.bins=unlist(lapply(my.vList,viewMins) )# max of each view
                
              }else{stop("wrong 'bin.op' option given")}
              
              # make a matrix from those views            
              mat=matrix( (sum.bins-target@add)/target@multiply, ncol=bin.num,byrow=TRUE)
              
              return(mat)
                              
            }else{
              
              if(bin.op=="mean")
              {
                sum.bins=unlist(lapply(my.vList,viewMeans) )# sum of each view
              }
              else if (bin.op=="max")
              {
                sum.bins=unlist(lapply(my.vList,viewMaxs) )# max of each view
              }
              else if (bin.op=="min")
              {
                sum.bins=unlist(lapply(my.vList,viewMins) )# max of each view
                
              }else{stop("wrong 'bin.op' option given")}
              
              # make a matrix from those views            
              mat=matrix( sum.bins, ncol=bin.num,byrow=TRUE)
              
              mat
              
            }      
 
}


#######################################
# S4 functions
#######################################



#' get scores overlapping with windows in a scoreMatrix object
#'
#' A scoreMatrix object can be used to draw average profiles or heatmap of read coverage or wig track-like data.
#' \code{windows} can be a predefined region such as CpG islands or gene bodies that are not necessarily equi-width.
#' Each window will be chopped to equal number of bins based on \code{bin.num} option.
#'
#' @param target a \code{RleList} or a \code{modRleList} object to be overlapped with ranges in \code{windows}
#' @param windows a \code{GRanges} object that will be randomly placed across the genome and overlap of these random regions with \code{target} will be the background distribution of association between \code{target} and \code{query}.
#' @param bin.num A single \code{integer} value denoting how many bins there should be for each window
#' @param bin.op A bin operation that is either one of the following strings: "max","min","mean". The operation is applied on the values in the bin. Defaults to "mean"
#' @param strand.aware If TRUE (default: FALSE), the strands of the windows will be taken into account in the resulting \code{scoreMatrix}. If the strand of a window is -, the values of the bins for that window will be reversed
#' @param ... parameters to be passed to \code{modCoverage} function
#'
#' @usage scoreMatrix(target,windows)
#' @return returns a \code{scoreMatrix} object
#' @export
#' @docType methods
#' @rdname scoreMatrixBin-methods           
setGeneric("scoreMatrixBin",function(target,windows,bin.num=10,bin.op="mean",strand.aware=FALSE,...) standardGeneric("scoreMatrixBin") )

#' @aliases scoreMatrixBin,GRanges,RleList-method
#' @rdname scoreMatrixBin-methods
setMethod("scoreMatrixBin",signature("RleList","GRanges"),
          function(target,windows,bin.num,bin.op,strand.aware){

            #check if windows lengths exceeds the length of feature based chromosomes
            r.chr.len=lapply(target,length)
            constraint=GRanges(seqnames=names(r.chr.len),IRanges(start=rep(1,length(r.chr.len)),end=unlist(r.chr.len))  )
            windows=subsetByOverlaps(windows, constraint,type = "within",ignore.strand = TRUE)
            
            if(! strand.aware)
            {
              mat=make.scoreMatrixBin(target, windows, bin.num, bin.op)
              new("scoreMatrix",mat)
            }else if(unique(strand(windows)) %in% "-" ){
              
              mat1=make.scoreMatrixBin(target,windows[strand(windows)=="-",],bin.num,bin.op)
              mat2=make.scoreMatrixBin(target,windows[strand(windows) != "-",],bin.num,bin.op)
              
              new("scoreMatrix",rbind( mat1[ncol(mat1):1],mat2) )
                
            }else{
              mat=make.scoreMatrixBin(target,windows,bin.num,bin.op)
              new("scoreMatrix",mat)              
            }             

      
})

#' @aliases  scoreMatrixBin,GRanges,GRanges,ANY-method
#' @rdname scoreMatrixBin-methods
setMethod("scoreMatrixBin",signature("GRanges","GRanges"),
          function(target,windows,bin.num,bin.op,strand.aware,...){
            
            
            #make coverage vector (modRleList) from target
            target.rle=modCoverage(target,...)
            
            # call scoreMatrix function
            scoreMatrixBin(target.rle,windows,bin.num,bin.op,strand.aware)
            
})
