#######################################
# S3 functions
#######################################
# Combine a scoreMatrix into a scoreMatrixList object

# ---------------------------------------------------------------------------- #
#  when a ScoreMatrixList is a first argument

c.ScoreMatrixList<-function(..., recursive = FALSE, use.names = TRUE) {
  y <- list(...)
  #Initialization of a ScoreMatrixList object
  p <- y[[1]]
  n <- length(p) 
  for (j in 2:length(y)){
    if(class(y[[j]]) == "ScoreMatrix"){
      #combine j-th scoreMatrix into a scoreMatrixList object
      p[[n+1]] <- y[[j]]
      ##add the labels of additional scoreMatrix objects
      #in case of lack of the label name, add a numbered "scoreMatrix" name
      if(identical(c(""), names(y)[[j]]) || is.null(names(y)[[j]])){
        names(p)[n+1] <- paste("scoreMatrix", n+1, sep = "")
        #add a label name
      }else{
        names(p)[n+1] <- names(y)[[j]]
      }
      n <- n + 1
    }else if(class(y[[j]]) == c("ScoreMatrixList")){
      for (i in 1:length(y[[j]])){
        #combine i-th scoreMatrix obtained from j-th scoreMatrixList into a scoreMatrixList object
        p[[n+1]] <- y[[j]][[i]]
        ##add the labels of additional scoreMatrix objects
        #in case of lack of the label name, add a numbered "scoreMatrix" name
        if(identical(c(""), names(y[[j]])[i]) || is.null(names(y[[j]])[i])){
          names(p)[n+1] <- paste("scoreMatrix", n+1, sep = "")
          #add a label name
        }else{
          names(p)[n+1] <- names(y[[j]])[i]
        }
        n <- n + 1
      }
    }
  }
  validObject(p)
  return(p)
}

# ---------------------------------------------------------------------------- #
#  when a ScoreMatrix is a first argument

c.ScoreMatrix<-function(..., recursive = FALSE, use.names = TRUE) {
  y <- list(...)
  n <- 1
  #Initialization of a ScoreMatrixList object
  p <- as(list(),"ScoreMatrixList")
  for (j in 1:length(y)){
    if(class(y[[j]]) == "ScoreMatrix"){
      #combine j-th scoreMatrix into a scoreMatrixList object
      p[[n]] <- y[[j]]
      ##add the labels of additional scoreMatrix objects
      #in case of lack of the label name, add a numbered "scoreMatrix" name
      if(identical(c(""), names(y)[[j]]) || is.null(names(y)[[j]])){
        names(p)[n] <- paste("scoreMatrix", n, sep = "")
        #add a label name
      }else{
        names(p)[n] <- names(y)[[j]]
      }
      n <- n + 1
    }else if(class(y[[j]]) == c("ScoreMatrixList")){
      for (i in 1:length(y[[j]])){
        #combine i-th scoreMatrix obtained from j-th scoreMatrixList into a scoreMatrixList object
        p[[n]] <- y[[j]][[i]]
        ##add the labels of additional scoreMatrix objects
        #in case of lack of the label name, add a numbered "scoreMatrix" name
        if(identical(c(""), names(y[[j]])[i]) || is.null(names(y[[j]])[i])){
          names(p)[n] <- paste("scoreMatrix", n, sep = "")
          #add a label name
        }else{
          names(p)[n] <- names(y[[j]])[i]
        }
        n <- n + 1
      }
    }
  }
  validObject(p)
  return(p)
}