#######################################
# S3 functions
#######################################
# Combine a scoreMatrix into a scoreMatrixList object

# ---------------------------------------------------------------------------- #
#  when a ScoreMatrixList is a first argument
c.ScoreMatrixList<-function(..., recursive = FALSE, use.names = TRUE) {
  y <- list(...)
  n <- length(y)
  y3 <- list()
  Ln <- vector()
  
  for (i in 1:n){
    #combine the scoreMatrixList object
    if (is(y[[i]], "ScoreMatrixList")){
      y2 <- as(y[[i]], "list")
      y3 <- append(y3, y2)
      #combine its label name
      Ln <- append(Ln, names(y[[i]]))
    #combine the scoreMatrix object
    }else if(is(y[[i]], "ScoreMatrix")){
      y3 <- append(y3, list(y[[i]]))
      #combine its label name
      if (identical(c(""), names(y[i])) || is.null(names(y[i]))) {
        Ln <- append(Ln, c(""))
      }else{
        Ln <- append(Ln, names(y[i]))
      } 
      
    }
  }  
  #create the ScoreMatrixList
  SML <- as(y3, "ScoreMatrixList")
  names(SML) <- Ln
  
  validObject(SML)
  return(SML)
}

# ---------------------------------------------------------------------------- #
#  when a ScoreMatrix is a first argument

c.ScoreMatrix<-function(..., recursive = FALSE, use.names = TRUE) {
  y <- list(...)
  n <- length(y)
  y3 <- list()
  Ln <- vector()
  
  for (i in 1:n){
    #combine the scoreMatrixList object
    if (is(y[[i]], "ScoreMatrixList")){
      y2 <- as(y[[i]], "list")
      y3 <- append(y3, y2)
      #combine its label name
      Ln <- append(Ln, names(y[[i]]))
   #combine the scoreMatrix object
    }else if(is(y[[i]], "ScoreMatrix")){
      y3 <- append(y3, list(y[[i]]))
      #combine its label name
      if (identical(c(""), names(y[i])) || is.null(names(y[i]))) {
        Ln <- append(Ln, c(""))
      }else{
        Ln <- append(Ln, names(y[i]))
      } 
      
    }
  }  
  #create the ScoreMatrixList
  SML <- as(y3, "ScoreMatrixList")
  names(SML) <- Ln
  
  validObject(SML)
  return(SML)
}

