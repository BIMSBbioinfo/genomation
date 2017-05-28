#######################################
# S3 functions
#######################################
# Combine a scoreMatrix into a scoreMatrixList object

# ---------------------------------------------------------------------------- #
#  when a ScoreMatrixList is a first argument
c.ScoreMatrixList<-function(..., recursive = FALSE, use.names = TRUE) {
  y <- list(...)
  
  #find indexes of both ScoreMatrix and ScoreMatrixList objects
  m <- sapply(y,function(x) is(x,"ScoreMatrix"))
  l <- sapply(y,function(x) is(x,"ScoreMatrixList"))
  m2 <- c(1:length(y))[m]
  l2 <- c(1:length(y))[l]
  
  #combine all ScoreMatrixList objects
  y2 <- y[l2]
  n <- length(y2)
  y3 <- list()
  Ln <- vector()
  for (i in 1:n){
    #combine the scoreMatrixList objects
    y2L <- as(y2[[i]], "list")
    y3 <- append(y3, y2L)
    #combine their names
    Ln <- append(Ln, names(y2[[i]]))
  }
  SML <- as(y3, "ScoreMatrixList")
  names(SML) <- Ln
  
  #when only ScoreMatrixList objects are added by a user (lack of the scoreMatrix object) --> stop here
  if (isEmpty(m2)){
    validObject(SML) 
    return(SML)
  }
  
  #combine all ScoreMatrix objects into the existing ScoreMatrixList
  p <- ScoreMatrixList(append(SML@.Data, y[m2]))
  
  #combine their names 
  y4 <- y[m2]
  n2 <- length(y4)
  Mn <- vector()
  for (j in 1:n2){
    if (identical(c(""), names(y4[j])) || is.null(names(y4[j]))) {
       Mn = append(Mn, c(""))
    }else{
      Mn = append(Mn, names(y4[j]))
    } 
  }
  
  p@names <- append(SML@names, Mn)

  validObject(p)
  return(p)
}

# ---------------------------------------------------------------------------- #
#  when a ScoreMatrix is a first argument

c.ScoreMatrix<-function(..., recursive = FALSE, use.names = TRUE) {
  y <- list(...)
  
  #find indexes of both ScoreMatrix and ScoreMatrixList objects
  m <- sapply(y,function(x) is(x,"ScoreMatrix"))
  l <- sapply(y,function(x) is(x,"ScoreMatrixList"))
  m2<-c(1:length(y))[m]
  l2<-c(1:length(y))[l]
  
  #combine all ScoreMatrixList objects
  y2<-y[l2]
  n <- length(y2)
  y3 <- list()
  Ln <- vector()
  for (i in 1:n){
    #combine the scoreMatrixList objects
    y2L <- as(y2[[i]], "list")
    y3 <- append(y3, y2L)
    #combine their names
    Ln <- append(Ln, names(y2[[i]]))
  }
  SML <- as(y3, "ScoreMatrixList")
  names(SML) <- Ln
  
  #combine all ScoreMatrix objects into the existing ScoreMatrixList
  p <- ScoreMatrixList(append(SML@.Data, y[m2]))
  
  #combine their names 
  y4 <- y[m2]
  n2 <- length(y4)
  Mn <- vector()
  for (j in 1:n2){
    if (identical(c(""), names(y4[j])) || is.null(names(y4[j]))) {
      Mn = append(Mn, c(""))
    }else{
      Mn = append(Mn, names(y4[j]))
    } 
  }
  
  p@names <- append(SML@names, Mn)
  
  validObject(p)
  return(p)
}
