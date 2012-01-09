setGeneric("scoreMatrixList",function(...) standardGeneric("scoreMatrixList"))

# constructor
setMethod("scoreMatrixList", signature(...),
          function(...){
            l = as.list(...)
            return(new("scoreMatrixList"), ...)
})

# --------------------------------- #
# Validator
scoreMatrixList.Check = function(sml){
	errors = character()
	# checks whether all matrices are of class scoreMatrix
	if(!all(unlist(lapply(l, function(x)class(x) == 'scoreMatrix'))))
		errors = paste(errors, 'All elements for scoreMatrixList need to be of class scoreMatrix', sep='\n')

	# checks whether all matrices are numeric
	if(!all(unlist(lapply(l, function(x)all(is.integer(x) | is.numeric(x)))))))
		errors = paste(errors, 'Not all matrices are of type integer or numeric', sep='\n')
		
	if(length(errors) == 0) TRUE else errors 
}




# --------------------------------- #
# plot functions for score matrix list

# setGeneric("heamapProfile", function(mat.list, ...))

# setMethod("heatmapProfile", signature("scoreMatrixList"),
		  # function(mat.list, ...){
		  
			# par(cex=(0.1 +0.1*nsamp/1000) * nsamp/500, mar=c(2,2,2,2), oma=c(1,1,1,1), cex.axis=1.5, cex.main=2)
            # layout(matrix(1:3, ncol=3), widths=c(5,1,1))
                     
            # image(x=0:ncol(mat), y=0:nrow(mat), z = t(mat.ord), col=mat.cols, useRaster=T, main=name , xlab="Positon", ylab="Sample")
            # AddSep(mat, rowsep[-length(rowsep)], "black")
                        
            # image(t(matrix(as.numeric(fact[ord.fact]), ncol=1)), col=key.cols, axes=F)

            # plots the annotation for each group
            # t = tab/2
            # t[2:nfac] = t[1:(nfac-1)] +t[2:nfac]
            # plot.new()
            # plot.window(xlim=c(0,2), ylim=c(0,nrow(mat)))
            # text(1, cumsum(t)+2^(1:nfac), levels(fact), cex=1.5)
		  
			
		  
		  # }
# )

