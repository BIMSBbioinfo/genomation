# ---------------------------------------------------------------------------- #
#' Plot functions for score matrix list
#' 
#' Functions plot a ScoreMatrixList object as a panel of heatmaps
#'
#' @param mat.list a \code{ScoreMatrixList} object
#' @param mat.cols colors to be used for plotting
#' @param xmarks an integer number to lable the thick marks on the x axis of each heatmap.
#'  By default it takes the values of -ncol/2, 0, ncol/2
#' @param ymarks a vector of that will lable the thick marks on the y axis
#' @param y.at a numeric vector that will specify the positions of the thick marks on the y axis
#' @param xcex, ycex an integer number which controls the character expansion on x and y axis
#' @param cex.main an integer number which controls the character expansion of the plot label
#' @param mar a vector of length 5 which controls the size of the margins. 
#' The order is the following: below, left, up, right, spacing between consecutive plots
#' @param use.names whether to use the names of the ScoreMatrixList object to label each plot
#' @param main whether to use the names of the ScoreMatrixList object to label each plot
#' @param xlab, ylab name to be used for the x/y axis 
#' @param ... other options (obselete for now)

#' @usage heatmapProfile(mat.list, mat.cols=NULL, ...)

#' @examples
#'  l = lapply(seq(20, 40,5), function(x)as(matrix(rpois(1000, x), ncol=25), 'scoreMatrix'))
#'  
#'  
#' @export
#' @docType methods
#' @rdname heatmapProfile-methods
setGeneric("heatmapProfile", 
           function(mat.list, mat.cols=NULL, xmarks=NULL, 
                    ymarks=NULL, y.at=NULL, xcex=1.5, 
                    ycex=1.5, cex.main=3, mar=NULL, 
                    use.names=T, xlab=NULL, ylab=NULL, ...)
             standardGeneric("heatmapProfile"))

#' @aliases heatmapProfile,ScoreMatrixList-method
#' @rdname heatmapProfile-methods
setMethod("heatmapProfile", signature(mat.list="ScoreMatrixList"),
          function(mat.list, mat.cols, xmarks, 
                   ymarks, y.at, xcex, ycex, 
                   cex.main, mar, use.names, xlab, ylab, ...){
            
            dims = unlist(lapply(mat.list, nrow))
            if(!length(unique(dims)) == 1)
              stop('ScoreMatrixList does not contain matrices with the same number of rows')
            
            
            # default matrix colors
            if(is.null(mat.cols))
              mat.cols = colorRampPalette(c('lightgray','darkblue'), 
                                          interpolate='spline')(20)
            
            
            # checks the margin parameter
            if(is.null(mar)){
              mar = rep(3, 5)
            }else if(! length(mar) == 5 | !all(is.numeric(mar))){
              stop('mar is not of length 5')
            }
            
            if(use.names)
              main = names(mat.list)
            
            
            # gets the dimension of the matrices
            ncols = unlist(lapply(mat.list, ncol))
            nrow = nrow(mat.list[[1]])
            len = length(mat.list)
            
            # sets the layout
            if(is.null(xlab) & is.null(ylab)){
              layout(matrix(1:len, ncol=len))
              
            }else if(! is.null(xlab) & is.null(ylab)){
              layout(matrix(c(2:(len+1), rep(1, len)), 
                            ncol=len, nrow=2, byrow=T), heights=c(20,1))
              .plotXYlab(xlab, 'x')
              
            }else if(is.null(xlab) & ! is.null(ylab)){
              layout(matrix(1:(len+1), ncol=len+1), widths=c(5,rep(20, len+1)))
              .plotXYlab(ylab, 'y')
              
            }else if(!is.null(xlab) & !is.null(ylab)){
              par(mar=c(0,1,0,1), oma=rep(0,4))
              layout(matrix(c(1,3:(len+2), 0,rep(2, len)), ncol=len+1, byrow=T), 
                     widths=c(4,rep(20, len+1)), 
                     heights=c(20,1))
              .plotXYlab(ylab, 'y')
              .plotXYlab(xlab, 'x')
            }
            
            # gets the tick marks labels and positions
            if(is.null(xmarks))
              xmarks = c(-ncols[1]/2, 0, ncols[1]/2)
            xpos = seq(1, ncols[1], length.out=length(xmarks))
            
            # sets the thick marks and labels on the y axis
            if(is.null(ymarks)){
              ymarks = round(fivenum(1:nrow))
            }
            if(is.null(y.at)){
              y.at = round(fivenum(1:nrow))
            }
            if(length(ymarks) != length(y.at))
              stop('ymarks and y.at do not have the same length')
            if(!is.null(y.at)){
              if(!is.numeric(y.at)){
                stop('y.at need to be a numeric variable')
              }
              if(any((y.at < 1) | (y.at > nrow))){
                stop('y.at values are outside of the matrix dimension')
              }
            }
            
            # plots the heatmaps
            for(i in 1:len){
              message('Plotting matrix: ', i,'\r')
              # sets the margins for each plot
              if(i == 1){	
                par(mar=c(mar[c(1,2,3)], mar[5]/2))
                
              }else if(i == len){
                par(mar=c(mar[1], mar[5]/2, mar[3], mar[4]))
              }else{
                par(mar=c(mar[1], mar[5]/2, mar[3], mar[5]/2))
              }
              
              image(x=1:ncols[1], y=1:nrow, z=t(mat.list[[i]]), 
                    main=main[i], cex.main=cex.main, col=mat.cols,
                    yaxt='n', xaxt='n', xlab='', ylab='', useRaster=T, ...)
              if(i==1){
                axis(2, at=y.at, las=2, cex.lab=2, labels=ymarks, cex.axis=ycex)
              }
              
              axis(1, at=xpos, labels=xmarks, cex.axis=xcex)
            }
            message('\nPlotting done\n')
            
          }
)





# ---------------------------------------------------------------------------- #
#' visual representation of ScoreMatrix using a heatmap 
#' The rows can be reordered using one factor and one numeric vector
#'
#' @param mat a \code{ScoreMatrix} object
#' @param fact a \code{factor} of length equal to \code{nrow(mat)}. Unused factor levels are dropped
#' @param ord.vec a \code{vector} of class \code{numeric} of the same length as mat, 
#'        which is going to be used for ordering of the rows
#' @param shift shift the start coordinate of the x axis (plot starts at -shift)
#' @param mat.cols a vector of colors used for plotting of the heatmap. 
#'        Default colors range from lightgray to darkblue.
#' @param fact.cols a vector of colors used for plotting of the factor key
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main plot name
#' @param class.names names for each factor class - has to have the same lenght as \code{levels(fact)}
#' @param ... other options to be passed to functions (obsolete at the moment)
#' @return nothing
#' 
#' @examples
#'   data(cage)
#'   data(promoters)
#'   myMat2=ScoreMatrix(target=cage,windows=promoters,
#'                         weight.col="tpm",strand.aware=TRUE)
#'   plot(colMeans(myMat2,na.rm=TRUE),type="l")
#'   heatMatrix(myMat2,fact=)
#' 
#' @docType methods
#' @rdname heatMatrix-methods
#' @export
setGeneric("heatMatrix", function(mat, 
                                  fact=NULL, 
                                  add.sep=TRUE, 
                                  ord.vec=NULL,
                                  shift=0, 
                                  mat.cols=NULL, 
                                  fact.cols=NULL, 
                                  xlab='Position', ylab='Region', 
                                  main='Positional profile', 
                                  class.names=NULL, use.names=FALSE, ...) 
  standardGeneric("heatMatrix") )

#' @aliases heatMatrix,ScoreMatrix-method
#' @rdname heatMatrix-methods
setMethod("heatMatrix", signature("ScoreMatrix"),
          function(mat, fact, add.sep, 
                   ord.vec, shift, mat.cols, 
                   fact.cols, xlab, ylab, 
                   main, class.names, use.names, ...){
            
            # -------------------------- #
            # parameter checking
            if(!is.null(fact) & length(fact) != nrow(mat))
              stop('Given factor does not have the same length as the matrix')
            if(!is.null(fact) & !is.factor(fact))
              stop('fact need to be an object of class factor')
            if(!is.null(ord.vec) & length(fact) != nrow(mat))
              stop('Given ordering vector does not have the same length as the matrix')
            if(!is.numeric(shift) | length(shift) > 1)
              stop('shift needs to be a numeric vector of length 1')
            # -------------------------- #
            # default values
            if(is.null(ord.vec))
              ord.vec = 1:nrow(mat)	
            
            if(is.null(mat.cols)){
              message('Using default mat.cols...\n')
              mat.cols = colorRampPalette(c('lightgray','darkblue'), 
                                          interpolate='spline')(20)
            }
            
            # -------------------------- #
            # plots the matrix
            AddSep = function(x, rowsep, col, sepwidth=c(0.05,0.5)){
              for(rsep in rowsep){
                rect(xleft =0, 
                     ybottom= (rsep), 
                     xright=ncol(x)+1,  
                     ytop = (rsep+1) - sepwidth[2], 
                     lty=1, lwd=1, col=col, border=col)
              }
            }
            
            if(!is.null(fact)){
              print('using factor')
              layout(matrix(c(1,2), ncol=2), widths=c(20,1))
              mat = mat[order(as.numeric(fact), ord.vec),]
            }else{
              mat = mat[ord.vec,]
            }
            
            par(mar=c(5,5,3,.1), oma=c(0,0,0,0))
            image(x=1:ncol(mat) - shift, y=1:nrow(mat), z=t(as.matrix(mat)), 
                  col=mat.cols, oma=c(0,0,0,0),
                  useRaster=TRUE, xlab=xlab, ylab=ylab, main=main, axes=FALSE)
            if(use.names==TRUE){
              axis(2, at=1:nrow(mat), labels=rownames(mat), las=2)
            }else{
              at = round(fivenum(1:nrow(mat)))
              axis(2, at=at, labels=at, las=2)
            }			
            
            if(!is.null(fact)){
              
              # drops unused levels from the factor
              fact = fact[1:length(fact), drop=TRUE]
              
              if(is.null(class.names))
                class.names = levels(fact)
              
              if(is.null(fact.cols)){
                message('Using default fact.cols...\n')
                fact.cols = getColors(length(levels(fact)))
              }
              
              classnum = table(fact)
              rowsep = cumsum(classnum)
              if(add.sep == TRUE)
                AddSep(mat, rowsep[-length(rowsep)], "black")
              
              par(mar=c(5,.5,3, max(nchar(class.names)/2)))
              image(x = 1:10, 
                    y = 1:nrow(mat), 
                    z=t(matrix(as.numeric(fact), nrow=length(fact), ncol=10)),
                    col = fact.cols, xaxt='n', yaxt='n', ylab='', xlab='',
                    oma=c(0,0,0,1))
              at = classnum/2
              at[-1] = at[-1] + at[-length(at)]
              at = cumsum(at)
              axis(side=4, at=at, labels=class.names, tick = F, las=2)
              
            }	
          }
)
