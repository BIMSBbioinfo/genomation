# functions for plotting metaProfiles:


# default color scheme for heatmaps

#.jets <- colorRampPalette(c("#334b8e", "#456ca7", "#71c6cd", "#8fc56c", 
#                           "#f3e92b","#e6762b","#d9272a"))
.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# winsorize a matrix using percentile ranges
.winsorize<-function(mat,rng){
  hi.th=quantile(mat,rng[2]/100,na.rm=TRUE)
  lo.th=quantile(mat,rng[1]/100,na.rm=TRUE)
  mat[mat>hi.th]=hi.th
  mat[mat<lo.th]=lo.th
  mat
}


# ---------------------------------------------------------------------------- #
#' Heatmap for meta-region profiles
#' 
#' Function calculates meta-profile(s) from a ScoreMatrix or a ScoreMatrixList, then
#' produces a heatmap or a set of stacked heatmaps for meta-region profiles
#'
#' @param mat \code{ScoreMatrix} or \code{ScoreMatrixList} to be plotted 
#' @param centralTend a character that determines central tendency of meta-profile(s). 
#'                     It takes "mean" (default) or "median".
#' @param profile.names a character vector for names of profiles. If NULL, 
#'                      the names
#'                      will be taken from names(mat) if mat is a 
#'                      \code{ScoreMatrixList} object.
#' @param xcoords a vector of numbers showing relative positions of the bases or 
#'                windows. It must match the number of columns in the \code{ScoreMatrix}
#'                For example: if there are 2001 elements in the matrices which
#'                are base-pair resolution 
#'                and they are centered around an anchor point like TSS, the xcoords
#'                argument should be -1000:1000. This argument is used to plot
#'                accurate x-axis labels for the plots.If NULL it will be equal
#'                to 1:ncol(mat).
#' @param meta.rescale if TRUE meta-region profiles are scaled to 0 to 1 range by
#'                     subracting the min from profiles and dividing them by 
#'                     max-min.
#' @param winsorize Numeric vector of two, defaults to c(0,100). This vector 
#'                  determines the upper and lower percentile values to limit the 
#'                  extreme values. For example, c(0,99) will limit the values to
#'                  only 99th percentile, everything above the 99 percentile will 
#'                  be equalized to the value of 99th percentile. This is useful 
#'                  for visualization of matrices that have outliers.             
#' @param col a vector of color pallete. 
#'        color scheme to be used. If NULL, a version of jet colors will be
#'            used.
#' @param legend.name a character label plotted next to the legend
#' @param cex.legend  A numerical value giving the amount by which 
#'                    legend axis marks should be magnified relative to the default
#' @param xlab label a character string for x-axis
#' @param main a character string for the plot title
#' @param cex.lab  A numerical value giving the amount by which 
#'                    axis labels (including 'legend.name') 
#'                    should be magnified relative to the default.
#' @param cex.axis  A numerical value giving the amount by which 
#'                    axis marks should be magnified relative to the default
#' 
#' @return returns meta-profile matrix invisibly.
#' 
#' 
#' @examples
#' # data(cage)
#' # data(promoters)
#' # scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE)
#' # data(cpgi)
#' # scores2=ScoreMatrix(target=cpgi,windows=promoters,strand.aware=TRUE)
#' 
#' # x=new("ScoreMatrixList",list(scores1,scores2))
#' # heatMeta(mat=x,legend.name="fg",cex.legend=0.8,main="fdf",cex.lab=6,
#' #          cex.axis=0.9)
#' @export
#' 
heatMeta<-function(mat, centralTend="mean",
                   profile.names=NULL,xcoords=NULL,col=NULL,
                   meta.rescale=FALSE, winsorize=c(0,100),
                   legend.name=NULL,cex.legend=1,xlab=NULL,
                   main="",cex.lab=1,cex.axis=1){
  
  # check class
  if(! class(mat) %in% c("ScoreMatrix","ScoreMatrixList"))
    stop("mat is not ScoreMatrix or ScoreMatrixList\n")
  # check centralTend
  if(! centralTend %in% c("median","mean"))
    stop("centralTend is not mean or median\n")
  
  
  # get meta profiles by taking the mean
  if( class(mat)=="ScoreMatrix" ){
    if(centralTend=="mean"){
      metas=list(colMeans(mat,na.rm=TRUE))
    }else{
      metas=list(apply(a, 2, function(x) median(x,na.rm=TRUE)))
    }
  }else if( class(mat)=="ScoreMatrixList" ){
    if(centralTend=="mean"){
    metas=lapply(mat,function(a) colMeans(a,na.rm=TRUE) )
    }else{
    metas=lapply(mat,function(a) apply(a, 2, function(x) median(x,na.rm=TRUE)) )
    }
  }  
  
  # if the ncols of matrices do not match do not plot anything
  if(length(unique(sapply(metas,length))) != 1){
    stop("ScoreMatrix number of columns do not match\n",
         "Try using binMatrix to make matrices with large number of columns",
         "equal to the smaller ones\n")
  }
  
  # get the default xcoordinates to plot
  if(!is.null(xcoords)){
    
    # if it is a two element vector of start and end coordinates 
    # of the first and the last column of the matrix
    if(length(xcoords)==2 & xcoords[1]<xcoords[2]){
      xcoords=seq(xcoords[1],xcoords[2],length.out=length(metas[[1]]) )
    }
    
    if(length(xcoords) != length(metas[[1]])) 
      stop("xcoords has wrong length: ",length(xcoords)," \n",
           "it should be equal to the number of columns of ScoreMatrices\n",
           "which is: ",length(metas[[1]]),"\n")    
  }else{
    xcoords=1:length(metas[[1]])
  }
  
  
  # try to get profile names from names of ScoreMatrixList
  if(is.null(profile.names) & !is.null(names(mat)) & class(mat)=="ScoreMatrixList" )
  {
    profile.names=names(mat)
  }
  # if user wants scaling
  if(meta.rescale){
    metas=lapply(metas,function(x) (x-min(x))/(max(x)-min(x)  )  )
  }
  
  img=do.call("rbind",metas)
  
  if(is.null(col)){
    col=.jets(100)
  }
  
  marOrg=par()$mar # get original parMar to be used later
  marNew=marOrg
  marNew[4]=6.1
  par(mar=marNew)
  image(x=xcoords,y=1:nrow(img),z=as.matrix(t(img[nrow(img):1,,drop=FALSE])),useRaster=TRUE,
        col=col,yaxt="n",ylab="",xlab=xlab,main=main,
        cex.lab=cex.lab,cex.axis=cex.axis)
  if(!is.null(profile.names)){
    axis(side=4,at=1:nrow(img),labels=rev(profile.names),
         las=2,cex.axis=cex.axis)
  }
  
  #plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure) # get the plot viewport from base graphics
  # showViewport(current.viewport());current.vpTree()
  #grid.text(c("one"),
  #          x=unit(1, "native"), y=unit(2, "native"),
  #          just="right", rot=60)
  
  if(par()$mar[2]<4.1){
    warning("left margin of the plot (set by mar in par()) should not be less then 4.1 lines")
  }
  # make view port for the legend
  legendVp <- viewport(width=unit(1, "lines"), height=unit(0.4, "npc"),
                       x = unit(3, "lines"), y = unit(0.5, "npc"),just="left")
  pushViewport(legendVp) # push the legend VP
  #current.viewport()
  current.vpTree()
  rng=range(img)
  .heatLegendY(min=rng[1],max=rng[2],cols=col,
               legend.name=legend.name,main=TRUE,cex.legend=cex.legend,
               cex.lab=cex.lab)
  popViewport(2) # remove the legend VP
  current.vpTree()         
  par(mar=marOrg)              
  
  invisible(metas)
}

#function based on plotrix::dispersion, 
#and additionally it takes into account NA values in data
#the problem is that when scoreMatrix is calculated then
#some of the columns have only one numerical value and rest of them is NA
#and then mean of such column is the value, but variation e.g. sd is NA,
.dispersion2 <- function(x,y,ulim,llim=ulim, border = NA, intervals=TRUE, ...){
  if (intervals) {
    llim <- y - llim
    ulim <- y + ulim
  }
  ulim.na <- is.na(ulim)
  if(any(ulim.na)){
    #selecting segments of numerical values that are located between NA's
    w <- which(ulim.na)
    previous.v <- 0 #location of NA to the left of segment
    for(i in 1:length(w)){ #for each segment
      if(w[i]==previous.v+1){
        #if there are some neighboring NA
        previous.v <- w[i]
        next
      }
      from <- previous.v+1
      to <- w[i]-1
      #running polygon separately for each segment
      polygon(c(x[from:to], rev(x[from:to])), 
              c(llim[from:to], rev(ulim[from:to])), 
              col = fill, border = border)      
      previous.v <- w[i]
    }
    from <- previous.v+1
    to <- length(x)
    polygon(c(x[from:to], rev(x[from:to])), 
            c(llim[from:to], rev(ulim[from:to])), 
            col = fill, border = border,
            lty=lty)
  }else{
    polygon(c(x, rev(x)), c(llim, rev(ulim)), col = fill, 
            border = border, ...)
  }
}

# ---------------------------------------------------------------------------- #
#' Line plot(s) for meta-region profiles
#' 
#' Function calculates meta-profile(s) from a ScoreMatrix or a ScoreMatrixList, then
#' produces a line plot or a set of line plots for meta-region profiles
#' 
#' @param mat \code{ScoreMatrix} or \code{ScoreMatrixList} object. If it is a 
#' \code{ScoreMatrixList} object, all matrices in the ScoreMatrixList should have 
#' the same number of 
#' columns.
#' @param centralTend a character that determines central tendency of meta-profile(s). 
#'                     It takes "mean" (default) or "median".
#' @param overlay If TRUE multiple profiles will be overlayed in the same plot
#'                (Default:TRUE). If FALSE, and mat is a ScoreMatrixList, consider
#'                using par(mfrow=c(1,length(mat)))  to see the plots from all
#'                matrices at once.
#' @param winsorize Numeric vector of two, defaults to c(0,100). This vector 
#'                  determines the upper and lower percentile values to limit the 
#'                  extreme values. For example, c(0,99) will limit the values to
#'                  only 99th percentile, everything above the 99 percentile will 
#'                  be equalized to the value of 99th percentile.This is useful 
#'                  for visualization of matrices that have outliers.
#' @param profile.names a character vector for names of the profiles. The order
#'        should be same as the as the order of ScoreMatrixList.
#' @param xcoords a numeric vector which designates 
#'        relative base positions of the meta-region profiles.
#'        For example, for a 2001 column ScoreMatrix, xcoord=-1000:1000 indicates
#'        relative positions of each column in the score matrix. If NULL (Default),
#'        xcoords equals to 1:ncol(mat) 
#' @param meta.rescale if TRUE meta-region profiles are scaled to 0 to 1 range by
#'                     subtracting the min from profiles and dividing them by max-min.
#'                     If dispersion is not FALSE, then dispersion will be scaled as well. 
#' @param smooth.func the function to smooth central tendency and dispersion bands (Default: NULL), e.g. 
#'                    stats::lowess. All NA's will be removed before smoothing.
#' @param line.col color of lines for \code{centralTend} of meta-region profiles. Defaults to colors from
#'        \code{rainbow()} function.
#' @param ylim same as \code{ylim} at \code{\link{plot}} function. 
#'             if NULL ylim is estimated from all meta-region profiles.
#' @param ylab same as \code{ylab} at \code{\link{plot}} function. 
#'             Default: "average score"
#' @param xlab same as \code{xlab} at \code{\link{plot}} function. 
#'             Default: "bases"
#' @param dispersion shows dispersion interval bands around \code{centralTend} (defualt:FALSE). It takes 
#'        one of the character:
#' \itemize{
#'  \item{"se"}{shows standard error of the mean and 95 percent confidence interval for the mean}
#'  \item{"sd"}{shows standard deviation and 2*(standard deviation)}
#'  \item{"IQR"}{shows 1st and 3rd quartile, and 
#'               confidence interval around the median based on the median +/- 1.57 * IQR/sqrt(n) (notches)}
#' }
#' @param dispersion.col color of band of \code{dispersion}.
#'        Defaults to colors from \code{rainbow()} function and transparency set to 0.5
#'        (rainbow(length(mat), alpha = 0.5)).
#' @param ... other options to \code{\link{plot}}
#' 
#' @return returns the meta-region profiles invisibly as a matrix.
#' 
#' @note
#' Notches show the 95 percent confidence interval for the median 
#' according to an approximation based on the normal distribution.
#' They are used to compare groups - if notches corresponding to adjacent base pairs
#' on the plot do not overlap, this is strong evidence that median differs.
#' Small sample sizes (5-10) can cause notches to extend beyond the interquartile range (IQR) 
#' (Martin Krzywinski \emph{et al}. \emph{Nature Methods 11}, 119-120 (2014))
#' 
#' TODO
#' If you'd like to visualize more than one ScoreMatrix (e.g. ScoreMatrixList containing two matrices)
#' Matrices are plotted in the same order as in the ScoreMatrixList. przydatna informacji przy
#' tym jak ktos zwizualizowac. Ustaw je w takiej kolejnosci w jakiej je chcesz zwizualizowac.
#' 
#' @examples
#' 
#' # data(cage)
#' # data(promoters)
#' # scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE)
#'
#' # data(cpgi)
#' # scores2=ScoreMatrix(target=cpgi,windows=promoters,strand.aware=TRUE)
#' 
#' # create a new ScoreMatrixList
#' # x=new("ScoreMatrixList",list(scores1,scores2))
#' # plotMeta(mat=x,overlay=TRUE,main="my plotowski")
#' 
#' # plot dispersion nd smooth central tendency and variation interval bands
#' # plotMeta(mat=x, centralTend="mean", dispersion="se", winsorize=c(0,99), 
#' #         main="Dispersion as interquartile band", lwd=4, 
#' #         smooth.func=function(x) stats::lowess(x, f = 1/5))
#' 
#' @export
#' @docType methods
#' @rdname plotMeta
#' 
plotMeta<-function(mat, centralTend="mean",
                   overlay=TRUE,winsorize=c(0,100),
                   profile.names=NULL,xcoords=NULL,
                   meta.rescale=FALSE,
                   smooth.func=NULL,
                   line.col=NULL,
                   dispersion=FALSE,dispersion.col=NULL,
                   ylim=NULL,ylab="average score",xlab="bases", ...){
  
  # check class
  if(! class(mat) %in% c("ScoreMatrix","ScoreMatrixList"))
    stop("mat is not ScoreMatrix or ScoreMatrixList\n")
  # check centralTend args
  if(! centralTend %in% c("median","mean"))
    stop("centralTend is not mean or median\n")
  # check dispersion args
  disp.args <- c("se","sd","IQR") #dispersion arguments
  if(! dispersion %in% c(disp.args,FALSE))
    stop("dispersion is not FALSE, 'se', 'sd' or 'IQR'\n")
  
  
  if(is.null(line.col) & dispersion==FALSE)
    line.col=ifelse(is.list(mat),
                    list(rainbow(length(mat))),
                    "black")[[1]]
  if(is.null(line.col) & dispersion!=FALSE & is.null(dispersion.col)){
    dispersion.col=ifelse(is.list(mat),
                          list(rainbow(length(mat), alpha = 0.4)),
                          rainbow(1, alpha=0.4))[[1]]
    line.col=ifelse(is.list(mat),
                    list(rainbow(length(mat))),
                    rainbow(1))[[1]]
  }
  
  #mat is always a list/ScoreMatrixList
  if( class(mat)=="ScoreMatrix" ){
    mat <- list(mat)
  }
  
  # if the ncols of matrices do not match do not plot anything
  if(length(unique(sapply(mat,length))) != 1){
    stop("ScoreMatrix number of columns do not match\n",
         "Try using binMatrix to make matrices with high number of columns",
         "equal\n")
  }
  
  #init of some variables before for loop
  if(dispersion %in% disp.args){
    bound2<-list()
    if(dispersion=="IQR"){
      q1<-list(); q3<-list();
    }else{
      bound1 <- list()
  }}
  metas<-list()
  
  
  for(i in 1:length(mat)){
    
    # this can set extreme values to given percentile
    if(winsorize[2]<100 | winsorize[1]>0){
      mat[[i]]=.winsorize(mat[[i]],winsorize)
    }
    
    # get meta profiles by taking the mean/median
    if(centralTend=="mean"){
      if(dispersion=="IQR"){
        warning("dispersion is set to show 1st and 3rd quartile and 
                confidence interval around the median, 
                but centralTend is 'mean'. Setting centralTend to 'median'..\n")
        metas[[i]]=colMedians(mat[[i]], na.rm=TRUE)
      }else{
        metas[[i]]=colMeans(mat[[i]], na.rm=TRUE) 
      }
    }else if(centralTend=="median"){
      if(dispersion=="se"){
        warning("dispersion is set to standard error of the mean and 95% confidence interval for the mean, but
                centralTend is 'median'. Setting centralTend to 'mean'\n")
        metas[[i]]=colMeans(mat[[i]],na.rm=TRUE) 
      }else{
        metas[[i]]=colMedians(mat[[i]], na.rm=TRUE)
      }
    }
    
    # calculate dispersion around the mean/median
    if(dispersion %in% disp.args){      
      if(dispersion=="se"){
        bound1[[i]] <- std.error(mat[[i]], na.rm = TRUE)
        bound2[[i]] <- bound1[[i]] * 1.96
      }else if(dispersion=="sd"){
        bound1[[i]] <- colSds(mat[[i]], na.rm=TRUE)
        bound2[[i]] <- bound1[[i]] * 2
      }else if(dispersion=="IQR"){
        q <- colQuantiles(mat[[i]], probs=c(0.25, 0.75), na.rm=TRUE)
        q1[[i]] <- q[,1] #1st quartile
        q3[[i]] <- q[,2] #3rd quartile
        n<-ncol(mat[[i]])
        bound2[[i]] <- (1.57*(q3[[i]] - q1[[i]])) / sqrt(n) #notch
      }
    }
    
    if(meta.rescale){
      val2unit <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))} 
      metas[[i]]=val2unit(metas[[i]])
      if(dispersion %in% disp.args){
        if(dispersion=="IQR"){
          bound1[[i]]=val2unit(bound1[[i]])
          bound2[[i]]=val2unit(bound2[[i]])
        }else{
          q1[[i]]=val2unit(q1[[i]])
          q3[[i]]=val2unit(q3[[i]]) 
          bound2[[i]]=val2unit(bound2[[i]])
        }
      }
    }
    
    #smoothing using function smooth.func defined by user
    if(!is.null(smooth.func)){
      #first removing NA's and then smoothing
      metas[[i]] <- smooth.func(metas[[i]])$y
      if(dispersion %in% disp.args){
        bound2[[i]] <- smooth.func(bound2[[i]])$y
        if(dispersion=="IQR"){
          q1[[i]] <- smooth.func(q1[[i]])$y
          q3[[i]] <- smooth.func(q3[[i]])$y
        }else{
          bound1[[i]] <- smooth.func(bound1[[i]])$y
        }
      }
    }
    
  }

  # get the default xcoordinates to plot
  if(!is.null(xcoords)){
    
    # if it is a two element vector of start and end coordinates 
    # of the first and the last column of the matrix
    if(length(xcoords)==2 & xcoords[1]<xcoords[2]){
      xcoords=seq(xcoords[1],xcoords[2],length.out=length(metas[[1]]) )
    }
    
    if(length(xcoords) != length(metas[[1]])) 
      stop("xcoords has wrong length: ",length(xcoords)," \n",
           "it should be equal to the number of columns of ScoreMatrices\n",
           "which is: ",length(metas[[1]]),"\n")    
  }else{
    xcoords=1:length(metas[[1]])
  }
  
  
  # if ylim is not NULL, change the ranges to plot to ylim
  if(!is.null(ylim)){
    myrange=ylim
  }else{
    myrange=range(unlist(metas),na.rm = TRUE)
    if(dispersion %in% disp.args){
      bound2.max <- max(unlist(bound2), na.rm = TRUE)
      myrange[2] <- myrange[2] + abs(bound2.max)
      myrange[1] <- myrange[1] - abs(bound2.max)
    }
  }
  
  marOrg=par()$mar # get original parMar to be used later
  marNew=marOrg
  marNew[4]=8
  par(mar=marNew) # extend right margin for the legend
  par(xpd=TRUE) # do this so that you can plot legend out of the plotting box
  
  if(overlay & length(metas)>1){
    # plot overlayed lines
    if(dispersion %in% disp.args){
      plot(xcoords,metas[[1]],type="l",col=dispersion.col[1],
           ylim=myrange,ylab=ylab,xlab=xlab,...)
      for(i in 1:length(metas) ){
        .dispersion2(xcoords, metas[[i]], bound2[[i]],
                   fill=dispersion.col[i],lty=1)
        if(dispersion=="IQR"){
          .dispersion2(xcoords, metas[[i]], llim=metas[[i]]-q1[[i]], ulim=q3[[i]]-metas[[i]], 
                       fill=dispersion.col[i],lty=1)
        }else{
          .dispersion2(xcoords, metas[[i]], bound1[[i]],
                     fill=dispersion.col[i],lty=1)
        }  
      }
      for(j in 1:length(metas) ){
        i <- d$index[j]
        lines(xcoords,metas[[i]],col=line.col[j],...)
      }
      
    }else{ #without plotting dispersion
      plot(xcoords,metas[[d$index[1]]],type="l",col=line.col[1],
           ylim=myrange,ylab=ylab,xlab=xlab,...)
      for(i in 2:length(metas) ){
        lines(xcoords,metas[[i]],col=line.col[i],...)
      }
    }
    
    # if profile names are given, plot them as legend
    if(!is.null(profile.names))
      if(dispersion %in% disp.args){
        legend(max(xcoords)+0.05*max(xcoords),myrange[2],legend=profile.names
               ,fill=dispersion.col[d$index],bty="n", border=line.col)
      }else{
        legend(max(xcoords)+0.05*max(xcoords),myrange[2],legend=profile.names
               ,fill=line.col,bty="n")
      }
  }else{ # plot things one by one, in this case user must use par
    
    for(j in 1:length(metas)){
      i <- d$index[j]
      plot(xcoords,metas[[i]],type="l",col=line.col[i],
           ylim=myrange,ylab=ylab,xlab=xlab,...)
      if(dispersion %in% disp.args){
        .dispersion2(xcoords, metas[[i]], bound2[[i]],col=dispersion.col[j],
                   fill=dispersion.col[j],lty=1)
        if(dispersion=="IQR"){
          .dispersion2(xcoords, metas[[i]], llim=metas[[i]]-q1[[i]], ulim=q3[[i]]-metas[[i]],
                     fill=dispersion.col[j],lty=1) 
        }else{
          .dispersion2(xcoords, metas[[i]], bound1[[i]],col=dispersion.col[j],
                     fill=dispersion.col[j],lty=1)
        }
      }
      lines(xcoords,metas[[i]],col=line.col[j],...)
    }
  }
  # revert par shit to its original state
  par(xpd=FALSE)
  par(mar=marOrg)
  
  invisible(do.call("rbind",metas))
}


# put a y axis legend
.heatLegendY<-function(min,max,cols,legend.name,main=TRUE,cex.legend=1,
                       cex.lab=1){
  
  # get value range as a vector of 100
  vals=seq(min,max,length.out=100)
  rng <- range(vals, na.rm = TRUE) # get min/mqx
  m <- (vals - min)/(max-min) # get normalized range
  rasta= rgb(colorRamp(cols)(m), maxColorValue= 255) # get color for each element of range
  
  grid.raster( rev(rasta), interpolate=FALSE,height = unit(1, "npc"),
               width=unit(1, "npc")) # make the legend
  # make legend ticks
  at = seq(0,1,length.out=5); label = seq(min,max,length.out=5)
  
  #make the axis of the legend
  grid.yaxis(at=at,label=formatC(label,digits=2,format="g"),main=main,
             edits=gEdit("labels", rot=90,hjust=0.5),
             gp=gpar(cex=cex.legend)) # make axis for legend
  my.x = -2
  if(main==FALSE) 
    my.x=3.4
  grid.text(legend.name,rot=90,x=unit(my.x, "npc"),gp=gpar(cex=cex.lab))
  
}

.heatLegendX<-function(min,max,cols,legend.name,main=TRUE,cex.legend=1,
                       cex.lab=1,hjust=0,vjust=0){
  
  # get value range as a vector of 100
  vals=seq(min,max,length.out=100)
  rng <- range(vals, na.rm = TRUE) # get min/mqx
  m <- (vals - min)/(max-min) # get normalized range
  rasta= rgb(colorRamp(cols)(m), maxColorValue = 255) # get color for each element of range
  
  grid.raster( matrix((rasta),nrow=1), interpolate=FALSE,height = unit(1, "npc"),
               width=unit(1, "npc")) # make the legend
  # make legend ticks
  label = pretty(c(min,max),n=5);at = seq(0,1,length.out=length(label)); 
  
  #make the axis of the legend
  grid.xaxis(at=at,label=label,main=main,
             edits=gEdit("labels",hjust=hjust,vjust=vjust),
             gp=gpar(cex=cex.legend)) # make axis for legend
  my.y = -3
  grid.text(legend.name,y=unit(my.y, "npc"),gp=gpar(cex=cex.lab))
  
}

# convert a matrix or vector to color matrix
# to be used in grid.raster()
.convertToColors <- function(mat,cols,rng=NULL) {
  
  if(is.null(rng)){
    # Produce 'normalized' version of matrix, with values ranging from 0 to 1
    rng <- range(mat, na.rm = TRUE)
  }
  
  m <- (mat - rng[1])/(diff(rng))
  # Convert to a matrix of sRGB color strings
  #m2 <- m; class(m2) <- "character"
  m2<-matrix("transparent",ncol=ncol(m),nrow=nrow(m))
  m2[!is.na(m)] <- rgb(colorRamp(cols)(m[!is.na(m)]), maxColorValue = 255)
  #m2[is.na(m)] <- "transparent"
  return(m2)
}

# make a heatmap from a given matrix using grid.raster()
.gridHeat<-function(mat,col,xcoords,xlab,cex.lab,cex.axis,angle=0,
                    hjust=0,vjust=0,rng=NULL){
  
  mat2=.convertToColors(mat,col,rng)
  ras=grid.raster(mat2,interpolate = FALSE, width= unit(1, "npc"),
                  height=unit(1, "npc"))
  
  # make legend ticks
  at = seq(0,1,length.out=5); label = seq(min(xcoords),max(xcoords)
                                          ,length.out=5)
  
  ax=grid.xaxis(at=at,label=formatC(label,digits=4,format="g"),
                edits=gEdit("labels", rot=angle,hjust=hjust,vjust=vjust),
                gp=gpar(cex=cex.axis)) # make axis for legend
  
  grid.text(xlab,y=unit(-2.5, "lines"),gp=gpar(cex=cex.lab)) 
  #grid.draw(ax)
}


.rowSideCol<-function(group.vector,group.names=NULL,group.col=NULL,cex.lab=1){
  
  if( is.null(group.col) ){
    cols=rainbow(length(unique(group.vector)))
    img=cols[factor(group.vector,levels=unique(group.vector))]
  }else{
    img=group.col[factor(group.vector,levels=unique(group.vector))]    
  }
  grid.raster(img,interpolate = FALSE, width= unit(1, "npc"),
              height=unit(1, "npc"))
  
  # segment heights calculated from group.vector
  # will be used to put group names in the middle of the segment
  segh=as.vector(table(factor(group.vector,levels=unique(group.vector))))
  name.coord=1-((cumsum(segh)-(segh/2))/sum(segh)) # NPC coord
  
  if( is.null(group.names)){
    grid.text(unique(group.vector), y=unit(name.coord,"npc"), 
              x = unit(-0.5, "lines"),
              gp=gpar(cex=cex.lab),just="right")
  }else{
    grid.text(group.names, y=unit(name.coord,"npc"),
              x = unit(-0.5, "lines"),
              gp=gpar(cex=cex.lab),just="right")    
  }
  
}

# ---------------------------------------------------------------------------- #
#' Draw a heatmap of a given ScoreMatrix object
#' 
#' 
#' The function makes a heatmap out of given \code{ScoreMatrix} object. If desired
#' it can use clustering using k-means and plot cluster color codes as a sidebar. 
#' In addition, user can define groups of rows using 'group' argument.
#' 
#' @param mat a \code{ScoreMatrix} object
#' @param grid  if TRUE, grid graphics will be used. if FALSE, base graphics
#'               will be used on the top level, so users can use par(mfrow)
#'               or par(mfcol) prior to calling the function. Default:FALSE              
#' @param col   a vector of colors, such as the ones created by heat.colors(10).
#'              If NULL (which is default), jet color scheme (common in matlab
#'              plots) will be used.
#' @param xcoords a vector of numbers showing relative positions of the bases or 
#'                windows. It must match the number of columns in the \code{ScoreMatrix}. 
#'                Alternatively, it could be a numeric vector of two elements. Such
#'                as c(0,100) showing the relative start and end coordinates of the first
#'                and last column of the \code{ScoreMatrix} object.
#'
#' @param group a list of vectors of row numbers or a factor. This grouping is
#'              used for rowside colors of the heatmap. If it is a list,
#'              each element of the list must be a vector of row numbers. Names
#'              of the elements of the list will be used as names of groups. 
#'              If \code{group} is a factor
#'              , it's length must match the number of rows of the matrix, and 
#'              factor levels will be used as the names of the groups in the plot.
#
#'               
#' @param group.col a vector of color names to be used at the rowside colors if
#'                  \code{group} argument is given or \code{kmeans=TRUE}
#' @param order    Logical indicating if the rows should be ordered or not 
#'                 (Default:FALSE). If \code{order=TRUE} the matrix will be ordered
#'                 with rowSums(mat) values in descending order. If kmeans=TRUE
#'                 or \code{group} argument is provided, first the groups/clusters
#'                 will be ordered in descending order of sums of rows then, everything
#'                 within the clusters will be ordered by sums of rows.
#' @param winsorize Numeric vector of two, defaults to c(0,100). This vector 
#'                  determines the upper and lower percentile values to limit the 
#'                  extreme values. For example, c(0,99) will limit the values to
#'                  only 99th percentile, everything above the 99 percentile will 
#'                  be equalized to the value of 99th percentile.This is useful 
#'                  for visualization of matrices that have outliers.
#' @param kmeans    Logical indicating if kmeans clustering should be done on the 
#'                  rows or not (Default:FALSE).
#' @param k     Defaults to 3. It designates the number of clusters to be returned
#'              by kmeans clustering.
#' @param main a character string for the plot title
#' @param legend.name a character label plotted next to the legend
#' @param cex.legend  A numerical value giving the amount by which 
#'                    legend axis marks should be magnified relative to the default
#' @param xlab label a character string for x-axis of the heatmap
#' @param cex.lab  A numerical value giving the amount by which 
#'                    axis labels (including 'legend.name') 
#'                    should be magnified relative to the default.
#' @param cex.main  A numerical value giving the amount by which 
#'                    plot title should be magnified  
#' @param cex.axis  A numerical value giving the amount by which 
#'                    axis marks should be magnified relative to the default
#' @param newpage logical indicating if \code{grid.newpage()} function should be
#'                invoked if \code{grid=TRUE}.
#'                
#' @return
#'  returns kmeans clustering result invisibly, if kmeans=TRUE
#'                                    
#' @examples
#' # data(cage)
#' # data(promoters)
#' # scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE,
#' #                   weight.col="tpm")
#'
#'
#'
#' # heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
#' #           xcoords=-1000:1000,
#' #           cex.legend=0.8,main="CAGE clusters on promoters",cex.lab=1,
#' #           cex.axis=0.9,grid=FALSE)
#'
#' # set.seed(1000)
#' # heatMatrix(mat=scores1,legend.name="tpm",winsorize=c(0,99),xlab="region around TSS",
#' #         xcoords=-1000:1000,kmeans=TRUE,k=3,
#' #         cex.legend=0.8,main="CAGE clusters on promoters",cex.lab=1,
#' #         cex.axis=0.9,grid=FALSE)
#'           
#' 
#' 
#' @export
#' @rdname heatMatrix
heatMatrix<-function(mat,grid=FALSE,col=NULL,xcoords=NULL,
                     group=NULL,group.col=NULL,order=FALSE,
                     winsorize=c(0,100),
                     kmeans=FALSE,k=3,
                     main="",legend.name=NULL,cex.legend=1,xlab=NULL,cex.main=1,
                     cex.lab=1,cex.axis=1,newpage=TRUE
){
  
  if( class(mat) !="ScoreMatrix" ){stop("'mat' is not a ScoreMatrix object\n")}
  
  mat2=mat@.Data # get the matrix class, some operations are not transitive
  
  # if this is changed, a rowSide color map will be drawn
  # setting kmeans or group.list will populate this vector
  # and this function will check on its value later on
  # to decide to plot rowside colors or not
  group.vector=NULL
  group.names=NULL # will be plotted as annotation if filled in later
  
  
  # this can set extreme values to given percentile
  # better for visualization
  # alternative is to take log or sth
  # but that should be done prior to heatMatrix
  if(winsorize[2]<100 | winsorize[1]>0){
    hi.th=quantile(mat2,winsorize[2]/100,na.rm=TRUE)
    lo.th=quantile(mat2,winsorize[1]/100,na.rm=TRUE)
    mat2[mat2>hi.th]=hi.th
    mat2[mat2<lo.th]=lo.th
    
  }
  
  # do kmeans if requested
  if(kmeans){
    
    # impute values if there are NAs
    if(any(is.na(mat2)) ) {
      mat3=impute.knn(mat2 ,k = 10, 
                      rowmax = 0.5, colmax = 0.8, 
                      maxp = 1500)$data
      clu=kmeans(mat3,c=k)
    }
    else{
      # cluster 
      clu=kmeans(mat2,c=k)
    }
    
    # get group.vector centers, will be used at ordering later
    group.vector=clu$cluster
    kcenters=clu$centers
    
    # order things by clusters only
    mat2=mat2[order(group.vector),]
    group.vector=group.vector[order(group.vector)]
    
    # if user wants to order
    if(order){
      
      # replicate center value for each cluster
      g.factor=factor(group.vector,levels=unique(group.vector))
      cent.val=rowSums(kcenters,na.rm=TRUE)[g.factor]
      
      # order by centers,cluster id, and do ordering within clusters
      my.order=order(-cent.val,group.vector,-rowSums(mat2,na.rm=TRUE))
      
      # commence the new order: Novus Ordo Seclorum
      mat2=mat2[my.order,]
      group.vector=group.vector[my.order]
    }
    
    group.names=unique(group.vector) # to be used for the rowSide colors
  }
  
  
  # check conditions of group.list
  # group must not have duplicated numbers
  # warn if total number group elements is below nrow(mat)
  if(!is.null(group)  & !kmeans){
    
    # if group is a list of rowids, row numbers for original windows argument
    if(is.list(group)){
      # group must not have duplicated numbers
      win.numbs=(lapply(group, function(x) unique(x) ))
      win.vec=unlist(win.numbs)
      if(any(table(unlist(win.numbs))>1)){
        stop("'group' containing a list must not have duplicated numbers\n")
      }
      row.ids=rownames(mat2) # get row ids from the matrix
      group.vector=rep(0,length(row.ids)) # make a starting vector of zeros
      for(i in 1:length(win.numbs)){
        group.vector[row.ids %in% win.numbs[[i]] ]=i # make group vector
      }
      
      # if list has names, take it as group.names
      # group.names will be plotted with rowSide colors
      if(!is.null(names(group))){
        group.names=names(group)
      }
      
      # stop if total number group elements is below nrow(mat)
      if(all(group.vector==0)){
        stop("None of the elements in 'group' are a part of rownames(mat) \n")
      }
      
      # warn if total number group elements is below nrow(mat)
      if(any(group.vector==0)){
        warning("Number of elements in 'group' argument is less then nrow(mat) \n",
                "Dropping rows from 'mat' that are not contained in 'group'\n")
        mat2=mat2[group.vector>0,]
        group.vector=group.vector[group.vector>0]
      }       
      
      
    }
    else if(is.factor(group)){
      
      # check if the length is same as nrow(mat)
      if( length(group) != nrow(mat)){
        stop("'group' is a factor, and its length should be equal to nrow(mat)\n")
      }
      # arrange factor levels by the order of appeareance
      group=factor(as.character(group),levels=as.character(unique(group)))
      group.names=levels(group)
      levels(group)=1:length(levels(group))
      group.vector=as.numeric(group)
    }else{
      stop("'group' must be a factor or a list\n")
    }
    
    
    
    # order things by clusters only
    mat2=mat2[order(group.vector),]
    group.vector=group.vector[order(group.vector)]
    
    
    
    if(order){
      my.order=order(group.vector,-rowSums(mat2,na.rm=TRUE))
      
      # commence the new order: Novus Ordo Seclorum
      mat2=mat2[my.order,]
      group.vector=group.vector[my.order]
      names(group.vector)=rownames(mat2)
    }
    
  }else if(order & !kmeans ){ # if only ordering is needed no group or clustering
    order.vector       =rep(1,nrow(mat2))
    names(order.vector)=rownames(mat2)
    order.vector       = order.vector[order(-rowSums(mat2,na.rm=TRUE))]
    mat2               =mat2[order(-rowSums(mat2,na.rm=TRUE)),]
    
  }
  
  
  # THE PLOTTING STARTS HERE with ordered mat2
  if(!grid){
    plot.new()
    vps <- baseViewports()
    pushViewport(vps$figure) # get the plot viewport from base graphics
    
  }else{
    if(newpage)
      grid.newpage()
  }
  
  # get the default/given xcoordinates to plot
  if(!is.null(xcoords) & is.vector(xcoords)){
    
    # if it is a two element vector of start and end coordinates 
    # of the first and the last column of the matrix
    if(length(xcoords)==2 & xcoords[1]<xcoords[2]){
      xcoords=seq(xcoords[1],xcoords[2],length.out=ncol(mat2) )
    }
    
    if(length(xcoords) != ncol(mat2) ) 
      stop("xcoords has wrong length: ",length(xcoords)," \n",
           " it should be equal to the number of columns of ScoreMatrix\n",
           " which is",ncol(mat2),"\n")    
  }else{
    xcoords=1:ncol(mat2)
  }
  
  # get heatcolor scale
  if(is.null(col)){
    
    col=.jets(100)
  }
  
  
  
  
  # make legend viewport
  legendVp <- viewport(width=unit(0.7, "lines"), height=unit(0.4, "npc"),
                       x = unit(0.71, "npc"), y = unit(0.5, "npc"),just="left")
  pushViewport(legendVp)
  #grid.rect()
  rng=range(mat2,na.rm=TRUE)
  # make the legend 
  .heatLegendY(min=rng[1],max=rng[2],
               col,legend.name,main=FALSE,cex.legend,
               cex.lab)
  popViewport()
  
  
  
  # make heatmap viewport
  heatHeightNPC=0.7
  heatVp <- viewport(width=unit(0.5, "npc"), height=unit(heatHeightNPC, "npc"),
                     x = unit(0.2, "npc"), y = unit(0.5, "npc"),just="left")
  pushViewport(heatVp) # push the heatmap VP
  .gridHeat(mat2,col,xcoords,xlab,cex.lab,cex.axis,hjust=0.5) # make the heat  
  #upViewport(1) # up one level current.vpTree()
  popViewport()
  
  
  # make side colors
  if(!is.null(group.vector)){
    sideVp <- viewport(width=unit(0.05, "npc"), height=unit(heatHeightNPC, "npc"),
                       x = unit(0.145, "npc"), y = unit(0.5, "npc"),just="left")
    pushViewport(sideVp) # push the viewport
    
    grid.rect()
    
    # make the rowside colors and labels
    .rowSideCol(group.vector,group.names=group.names,
                group.col=group.col,cex.lab=cex.lab)
    
    popViewport()
  }
  
  # make the title
  #title.y=convertY(unit(0.85, "npc"), "lines")
  title.y=unit(0.9, "npc")
  grid.text(main, y=title.y, 
            x = unit(0.45, "npc"),
            gp=gpar(cex=cex.main))
  
  #return groups if k-means==TRUE
  if(!grid)popViewport()
  
  if(kmeans | !is.null(group.vector)){
    return(invisible(group.vector)) 
  }else if(order & is.null(group.vector) ){
    return(invisible(order.vector)) 
  }

}


# ---------------------------------------------------------------------------- #
#' Draw multiple heatmaps from a ScoreMatrixList object
#' 
#' The function plots multiple heatmaps for a ScoreMatrixList object side by side.
#' Each matrix can have different color schemes but it is essential that each matrix
#' is obtained from same regions or neighbouring regions. 
#' 
#' @param sml a \code{ScoreMatrixList} object
#' @param grid  if TRUE, grid graphics will be used. if FALSE, base graphics
#'               will be used on the top level, so users can use par(mfrow)
#'               or par(mfcol) prior to calling the function. Default:FALSE   
#' @param col    a color palette or list of color palettes, such as
#'               list(heat.colors(10),topo.colors(10)). If it is a list,
#'               it is length must match the number of matrices to be plotted.
#'               If it is a single palette
#'               every heatmap will have the same colors.
#' @param xcoords a vector of numbers showing relative positions of the bases or 
#'                windows or a list of vectors. 
#'                The elements of the list must match the number of columns in the
#'                corresponding \code{ScoreMatrix}. 
#'                Alternatively, the elements could be a numeric vector of two elements. Such
#'                as c(0,100) showing the relative start and end coordinates of the first
#'                and last column of the \code{ScoreMatrix} object. The remaining
#'                coordinates will be automatically matched in this case. If the
#'                argument is not a list but a single vector, then all heatmaps
#'                will have the same coordinate on their x-axis.
#'                
#' @param group a list of vectors of row numbers or a factor. The rows will be 
#'              reordered to match their grouping. The grouping is
#'              used for rowside colors of the heatmap. If it is a list,
#'              each element of the list must be a vector of row numbers. Names
#'              of the elements of the list will be used as names of groups. 
#'              If \code{group} is a factor
#'              , it's length must match the number of rows of the matrix, and 
#'              factor levels will be used as the names of the groups in the plot.
#
#'               
#' @param group.col a vector of color names to be used at the rowside colors if
#'                  \code{group} argument is given or \code{kmeans=TRUE}
#' @param order    Logical indicating if the rows should be ordered or not 
#'                 (Default:FALSE). If \code{order=TRUE} the matrix will be ordered
#'                 with rowSums of all matrices in descending order. If kmeans=TRUE
#'                 or \code{group} argument is provided, first the groups/clusters
#'                 will be ordered in descending order by the sums of rows, then everything
#'                 within the clusters will be ordered by the sums of rows.
#'
#' @param winsorize Numeric vector of two, defaults to c(0,100). This vector 
#'                  determines the upper and lower percentile values to limit the 
#'                  extreme values. For example, c(0,99) will limit the values to
#'                  only 99th percentile for a matrix, 
#'                  everything above the 99 percentile will 
#'                  be equalized to the value of 99th percentile.This is useful 
#'                  for visualization of matrices that have outliers.
#' @param kmeans    Logical indicating if kmeans clustering should be done on the 
#'                  rows or not (Default:FALSE).
#' @param k     Defaults to 3. It designates the number of clusters to be returned
#'              by kmeans clustering.
#' @param column.scale Logical indicating if matrices should be scaled or not,
#'                     prior to k-means clustering or ordering. Setting this
#'                     to TRUE scales the columns of the 
#'                     matrices using,
#'                     \code{scale()} function. scaled columns are only used for
#'                     clustering or ordering. Original scores are displayed for heatmaps.
#' @param matrix.main a vector of strings for the titles of the heatmaps. If NULL
#'                    titles will be obtained from names of the \code{ScoreMatrix}
#'                    objects in the \code{ScoreMatrixList} objects.
#' @param common.scale if TRUE (Default:FALSE) all the heatmap colors will be 
#'                    coming from the same score
#'                    scale, although each heatmap color scale can be different.
#'                    The color intensities will be coming from the same scale. 
#'                    The scale will be
#'                    determined by minimum of all matrices and maximum of all 
#'                    matrices. This is useful when all matrices are on the same
#'                    score scale. If FALSE, the color scale will be determined
#'                    by minimum and maximum of each matrix individually.
#' @param legend      if TRUE and color legend for the heatmap is drawn.                   
#' @param legend.name a vector of legend labels to be plotted with legends
#'                    of each heatmap. If it is a length 1 vector, all heatmaps
#'                    will have the same legend label.
#' @param cex.legend  A numerical value giving the amount by which 
#'                    legend axis marks should be magnified relative to the default
#' 
#' @param xlab  a vector of character strings for x-axis labels of the heatmaps. 
#'              if it is length 1, all heatmaps will have the same label.
#' @param cex.lab  A numerical value giving the amount by which 
#'                    axis labels (including 'legend.name') 
#'                    should be magnified relative to the default.
#' @param cex.main  A numerical value giving the amount by which 
#'                    plot title should be magnified  
#' @param cex.axis  A numerical value giving the amount by which 
#'                    axis marks should be magnified relative to the default
#' @param newpage logical indicating if \code{grid.newpage()} function should be
#'                invoked if \code{grid=TRUE}.
#' 
#' @return
#'  invisibly returns the order of rows, if kmeans=TRUE and/or order=TRUE
#'  
#' @examples
#' 
#' # data(cage)
#' # data(promoters)
#' # scores1=ScoreMatrix(target=cage,windows=promoters,strand.aware=TRUE)
#' 
#' # data(cpgi)
#' # scores2=ScoreMatrix(target=cpgi,windows=promoters,strand.aware=TRUE)
#' 
#' # sml=new("ScoreMatrixList",list(a=scores1,b=scores2))

#' # multiHeatMatrix(sml,kmeans=TRUE,k=2,matrix.main=c("cage","CpGi"),cex.axis=0.8)
#' 
#' # use with K-means
#' # multiHeatMatrix(sml,kmeans=TRUE,k=2,cex.axis=0.8,xcoords=c(-1000,1000),
#' #                 winsorize=c(0,99),
#' #                 legend.name=c("tpm","coverage"),xlab="region around TSS")
#' 
#' # use different colors
#' # require(RColorBrewer)
#' # col.cage= brewer.pal(9,"Blues")
#' # col.cpgi= brewer.pal(9,"YlGn")
#' # multiHeatMatrix(sml,kmeans=TRUE,k=2,cex.axis=0.8,xcoords=c(-1000,1000),
#' #                 winsorize=c(0,99),col=list(col.cage,col.cpgi),
#' #                 legend.name=c("tpm","coverage"),xlab="region around TSS")
#' 
#' 
#' 
#' 
#' @export
#'
multiHeatMatrix<-function(sml,grid=TRUE,col=NULL,xcoords=NULL,
                          group=NULL,group.col=NULL,order=FALSE,
                          winsorize=c(0,100),
                          kmeans=FALSE,k=3,column.scale=TRUE,
                          matrix.main=NULL,
                          common.scale=FALSE,legend=TRUE,
                          legend.name=NULL,cex.legend=0.8,
                          xlab=NULL,
                          cex.lab=1, cex.main=1,cex.axis=0.8,newpage=TRUE)
{
  
  if( class(sml) !="ScoreMatrixList" ){stop("'sml' is not a ScoreMatrix object\n")}
  
  # check if col and xcoords are lists
  # if so their number should be equal to the number of matrices
  if(is.list(col) & length(col) != length(sml)){
    stop("'col' is a list and its length does not match the length of ScoreMatrixList\n")
    col=NULL
  }
  
  if( is.list(xcoords) & length(xcoords) != length(sml) ){
    stop("'xcoords' is a list and its length does not match the length of ScoreMatrixList\n")
    xcoords=NULL
  }
  
  # check legend.name argument if , it is a vector
  if(!is.null(legend.name) & length(legend.name)>1 & length(legend.name) != length(sml)){
    stop("'legend.name' should match the length of the 'sml' ",
         "if it is not a length 1 vector\n")
  }else if (length(legend.name)==1){
    legend.name=rep(legend.name,length(sml))
  }
  
  # check xlab argument and set it
  if(!is.null(xlab) & length(xlab)>1 & length(xlab) != length(sml)){
    stop("'xlab' should match the length of the 'sml' ",
         "if it is not a length 1 vector\n")
  }else if (length(xlab)==1){
    xlab=rep(xlab,length(sml))
  } 
  
  #check number of rows if they match or not
  # if not, warn and run unionScoreMatrixList
  if( length(unique(sapply(sml,nrow)))>1  ){
    warning("\nThe row numbers are different\n",
            "attempting to get common rows and to reorder rows\n",
            "using 'intersectScoreMatrixList()'\n")
    sml=intersectScoreMatrixList(sml,reorder=TRUE)
  }
  
  
  
  mat.list=lapply(sml,function(x) x) # get the matrix class, some operations are not transitive
  
  
  # if this is changed, a rowSide color map will be drawn
  # setting kmeans or group.list will populate this vector
  # and this function will check on its value later on
  # to decide to plot rowside colors or not
  group.vector=NULL
  group.names=NULL # will be plotted as annotation if filled in later
  
  
  # this can set extreme values to given percentile
  # better for visualization
  # alternative is to take log or sth
  # but that should be done prior to heatMatrix
  if(winsorize[2]<100 | winsorize[1]>0){
    
    mat.list=lapply(mat.list,function(x) .winsorize(x,winsorize) )
    
  }
  
  
  # if order | kmeans is true
  # make a one large matrix by cbind
  if(kmeans | order){
    mat2=do.call("cbind",mat.list)
    if(column.scale){
      mat2=scale(mat2)
      mat2[is.nan(mat2)]=0
    }
  }
  
  # do kmeans if requested
  if(kmeans){
    
    
    # impute values if there are NAs
    if(any(is.na(mat2)) ) {
      mat3=impute.knn(mat2 ,k = 10, 
                      rowmax = 0.5, colmax = 0.8, 
                      maxp = 1500)$data
      clu=kmeans(mat3,c=k)
    }
    else{
      # cluster 
      clu=kmeans(mat2,c=k)
    }
    
    # get group.vector centers, will be used at ordering later
    group.vector=clu$cluster
    kcenters=clu$centers
    
    # order things by clusters only
    mat.list=lapply(mat.list,function(x) x[order(group.vector),])
    group.vector=group.vector[order(group.vector)]
    
    # if user wants to order
    if(order){
      
      # replicate center value for each cluster
      g.factor=factor(group.vector,levels=unique(group.vector))
      cent.val=rowSums(kcenters,na.rm=TRUE)[g.factor]
      
      # order by centers,cluster id, and do ordering within clusters
      my.order=order(-cent.val,group.vector,-rowSums(mat2,na.rm=TRUE))
      
      # commence the new order: Novus Ordo Seclorum
      mat.list=lapply(mat.list,function(x) x[my.order,])
      group.vector=group.vector[my.order]
    }
    
    group.names=unique(group.vector) # to be used for the rowSide colors
  }
  
  
  # check conditions of group.list
  # group must not have duplicated numbers
  # warn if total number group elements is below nrow(mat)
  if(!is.null(group)  & !kmeans){
    
    # if group is a list of rowids, row numbers for original windows argument
    if(is.list(group)){
      # group must not have duplicated numbers
      win.numbs=(lapply(group, function(x) unique(x) ))
      win.vec=unlist(win.numbs)
      if(any(table(unlist(win.numbs))>1)){
        stop("'group' containing a list must not have duplicated numbers\n")
      }
      row.ids=rownames(mat.list[[1]]) # get row ids from the matrix
      group.vector=rep(0,length(row.ids)) # make a starting vector of zeros
      for(i in 1:length(win.numbs)){
        group.vector[row.ids %in% win.numbs[[i]] ]=i # make group vector
      }
      
      # if list has names, take it as group.names
      # group.names will be plotted with rowSide colors
      if(!is.null(names(group))){
        group.names=names(group)
      }
      
      # stop if total number group elements is below nrow(mat)
      if(all(group.vector==0)){
        stop("None of the elements in 'group' are a part of rownames(mat) \n")
      }
      
      # warn if total number group elements is below nrow(mat)
      if(any(group.vector==0)){
        warning("Number of elements in 'group' argument is less then nrow(mat) \n",
                "Dropping rows from 'mat' that are not contained in 'group'\n")
        
        mat.list=lapply(mat.list,function(x) x[group.vector>0,])
        group.vector=group.vector[group.vector>0]
      }       
      
      
    }
    else if(is.factor(group)){
      
      # check if the length is same as nrow(mat)
      if(length(group) != nrow(mat.list[[1]])){
        stop("'group' is a factor, and its length should be equal to nrow(mat)\n")
      }
      # arrange factor levels by the order of appeareance
      group=factor(as.character(group),levels=as.character(unique(group)))
      group.names=levels(group)
      levels(group)=1:length(levels(group))
      group.vector=as.numeric(group)
    }else{
      stop("'group' must be a factor or a list\n")
    }
    
    
    
    # order things by clusters only
    mat.list=lapply(mat.list,function(x) x[order(group.vector),])
    group.vector=group.vector[order(group.vector)]
    
    
    
    if(order){
      
      # get cbound matrix for ordering
      mat2=do.call("cbind",mat.list)
      if(column.scale){
        mat2=scale(mat2)
        mat2[is.nan(mat2)]=0
      }
      my.order=order(group.vector,-rowSums(mat2,na.rm=TRUE))
      
      # commence the new order: Novus Ordo Seclorum
      mat.list=lapply(mat.list,function(x) x[my.order,])
      group.vector=group.vector[my.order]
    }
    
  }else if(order & !kmeans ){ # if only ordering is needed no group or clustering
    
    # get cbound matrix for ordering
    mat2=do.call("cbind",mat.list)
    if(column.scale){
      mat2=scale(mat2)
      mat2[is.nan(mat2)]=0
    }
    order.vector       =rep(1,nrow(mat2))
    names(order.vector)=rownames(mat2)
    my.order=order(-rowSums(mat2,na.rm=TRUE))
    mat.list=lapply(mat.list,function(x) x[my.order,])
    order.vector       = order.vector[my.order]
    
  }
  

  
  # THE PLOTTING STARTS HERE with ordered mat.list
  if(!grid){
    plot.new()
    vps <- baseViewports()
    pushViewport(vps$figure) # get the plot viewport from base graphics
    
  }else{
    if(newpage)
      grid.newpage()
  }
  
  # try to get matrix names from names of ScoreMatrixList
  if(is.null(matrix.main) & !is.null(names(sml)) & class(sml)=="ScoreMatrixList" )
  {
    matrix.main=  names(sml)
  }else if(!is.null(matrix.main) & length(matrix.main) != length(sml)){
    warning("'matrix.main' length does not match to the 'sml' length\n",
            "setting it to NULL")
    matrix.main=NULL
    
  }
  
  # calculate the width of heatmaps and everything else relative to that
  l.sml=length(mat.list) # get length of matrix
  # get the npc width of the heatmap
  # this will be our rudder when plotting
  # every coordinate and dimension will be relative to this value
  hw=40*(1-(convertX( unit(4,"lines"),"npc",valueOnly=TRUE)) )/(l.sml*44+7)
  hh=0.7 # this heatmap height in NPC
  hy=0.6 # heatmap y coordinate on the plot
  
  # if group.vector is not NULL, plot the rowside colors
  # make side colors
  if(!is.null(group.vector)){
    sideVp <- viewport(width=unit(hw/8, "npc"), height=unit(hh, "npc"),
                       x = convertX( unit(4,"lines"),"npc"), 
                       y = unit(hy, "npc"),just="left")
    pushViewport(sideVp) # push the viewport
    
    grid.rect()
    
    # make the rowside colors and labels
    .rowSideCol(group.vector,group.names=group.names,
                group.col=group.col,cex.lab=cex.lab)
    
    popViewport()
  }
  
  # heatmap start coordinates
  heat.startCoord=convertX( unit(4,"lines"),"npc",valueOnly=TRUE) + (hw/8) + (hw/20)
  
  # if the same scale to be used for all plots
  common.range=NULL
  if(common.scale){
    common.range=range(mat.list,na.rm=TRUE)
  }
  
  for(i in 1:length(mat.list)){
    
    
    # cxcoords stands for current xcoords
    if(is.list(xcoords)){
      cxcoords=xcoords[[i]]
    }else if( is.vector(xcoords) | is.null(xcoords) ){
      cxcoords=xcoords   
    }
    else if(!is.null(xcoords)){
      warning("xcoords should be a vector or a list or NULL,nothing else!!\n",
              "setting xcoords to NULL")
      cxcoords=NULL
    }
    
    # if it is a two element vector of start and end coordinates 
    # of the first and the last column of the matrix
    if(length(cxcoords)==2){
      cxcoords=seq(cxcoords[1],cxcoords[2],length.out=ncol(mat.list[[i]]))
    }
    
    # ccol: stands for current.color
    if(is.list(col)){
      ccol=col[[i]]
    }else if(is.vector(col) | is.null(col) ){
      ccol=col
    }else if(!is.null(col)){
      warning("col should be a vector or a list or NULL,nothing else!!\n",
              "setting colors to default")
      ccol=NULL
    }
    
    
    # get the default/given xcoordinates to plot
    if(!is.null(cxcoords)){
      if(length(cxcoords) != ncol( mat.list[[i]] ) ) 
        stop("xcoords has wrong length: ",length(cxcoords)," \n",
             " it should be equal to the number of columns of ScoreMatrix\n",
             " which is",ncol(mat.list[[i]]),"\n")    
    }else{
      cxcoords=1:ncol( mat.list[[i]] )
    }
    
    # get heatcolor scale
    if(is.null(ccol)) 
      ccol=.jets(100)
    
    
    # make heatmap viewport
    heatVp <- viewport(width=unit(hw, "npc"), height=unit(hh, "npc"),
                       x =unit(heat.startCoord,"npc"), 
                       y = unit(hy, "npc"),just="left")
    
    
    pushViewport(heatVp) # push the heatmap VP
    grid.rect()
    .gridHeat(mat.list[[i]],ccol,cxcoords,xlab[i],cex.lab,cex.axis,
              angle=60,hjust=0.6,vjust=-0.5,rng=common.range) # make the heat  
    #upViewport(1) # up one level current.vpTree()
    popViewport()
    
    if(legend){
      # make the legend viewport
      legendVp <- viewport(width=unit(hw*0.7, "npc"), 
                           height=unit(0.5, "lines"),
                           x =unit(heat.startCoord+hw*0.15,"npc"), 
                           y = unit(0.1, "npc"),
                           just="left")
      pushViewport(legendVp)
      grid.rect()
      
      if(common.scale){
        rng=common.range
      }else{  
        rng=range(mat.list[[i]],na.rm=TRUE)
      }
      # make the legend 
      .heatLegendX(min=rng[1],max=rng[2],
                   ccol,legend.name[i],main=TRUE,cex.legend,
                   cex.lab,vjust=-0.5,hjust=0.5)
      popViewport()
    }
    
    if(!is.null(matrix.main)){
      title.y=unit(0.96, "npc")
      grid.text(matrix.main[i], y=title.y, 
                x = unit(heat.startCoord+hw*0.5,"npc"),
                gp=gpar(cex=cex.main),just="bottom")
    }
    
    # update start coord for next one
    heat.startCoord=heat.startCoord+ hw+ (hw/20)
    
  }
  
  if(!grid){
    popViewport()
  }
  
  if(kmeans | !is.null(group.vector)){
    return(invisible(group.vector)) 
  }else if(order & is.null(group.vector) ){
    return(invisible(order.vector)) 
  }
  
}
