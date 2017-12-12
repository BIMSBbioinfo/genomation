#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp;

//#-------------------------------------------------------------------------#
//' Function that computes a mean value
//'
//' @param x NumericVector
//' @keywords internal
// [[Rcpp::export]]
double Mean_c(NumericVector x){
  //remove NA from the vector
  NumericVector x2 = na_omit(x);
  
  int sz2 = x2.size() ;// get the length of the input vector after removing NAs

  //return NA when the vector consists of NAs (after removing NAs the vector size is equal to 0)
  if(sz2 == 0){
    return NA_REAL;
  }
  
  return std::accumulate(x2.begin(), x2.end(), 0.0) / sz2 ; //calculate the mean value of the vector
}

//' Function that computes a mean value for each bin
//'
//' @param x NumericVector
//' @param n intiger - number of bins
//' @keywords internal
// [[Rcpp::export]]
NumericVector binMean(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  
  // if the bins number larger than vector size, return zeros 
  if(sz < n){
    if(all(is_na(x))){  // if the vector consists of NAs, return NAs
      res=rep(NA_REAL,n);
    }
    return res;
  }
  
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size,set the window size to 1
  if(sz == n){
    w_size = 1;
  }
  
  double prev = 0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i< n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    if(i == (n-1)){ // for the last bin
      end2 = sz;
    }
    
    prev = prev + w_size; // update the begining index of the slice
    
    NumericVector vec(&x[prev2], &x[end2]);
    res[i] = Mean_c(vec);
  }
  
  return res;
}

//#-------------------------------------------------------------------------#
//' Function that computes a median value
//'
//' @param x NumericVector
//' @keywords internal
// [[Rcpp::export]]
double Median_c(NumericVector x){
  //remove NA from the vector
  NumericVector x2 = na_omit(x);
  
  int sz2 = x2.size() ;// get the length of the input vector after removing NAs
  
  //return NA when the vector consists of NAs (after removing NAs the vector size is equal to 0)
  if(sz2 == 0){
    return NA_REAL;
  }
  
  if(sz2%2 == 0){
    std::sort(x2.begin(), x2.end());
    return (x2[(sz2/2)-1] + x2[sz2/2] ) / 2; 
  }else{
    std::nth_element(x2.begin(), x2.begin()+sz2/2, x2.end());
    return x2[sz2/2];
  }
  
}

//' Function that computes a median value for each bin
//'
//' @param x NumericVector
//' @param n intiger - number of bins
//' @keywords internal
// [[Rcpp::export]]
NumericVector binMedian(NumericVector x, int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  
  // if the bins number larger than vector size, return zeros 
  if(sz < n){
    if(all(is_na(x))){  // if the vector consists of NAs, return NAs
      res=rep(NA_REAL,n);
    }
    return res;
  }
  
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size,set the window size to 1
  if(sz == n){
    w_size = 1;
  }
  
  double prev=0; // index for start positions over vector
  int prev2 ;  // integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    if(i == (n-1)){ // for the last bin
      end2 = sz;
    }
    prev = prev + w_size; // update the begining index of the slice
    
    NumericVector vec(&x[prev2], &x[end2]);
    res[i] = Median_c(vec);
    
  }
  
  return res;
}



//#-------------------------------------------------------------------------#
//' Function that computes a max value
//'
//' @param x NumericVector
//' @keywords internal
// [[Rcpp::export]]
double Max_c(NumericVector x){
  //remove NA from the vector
  NumericVector x2 = na_omit(x);
  
  int sz2 = x2.size() ;// get the length of the input vector after removing NAs
  
  //return NA when the vector consists of NAs (after removing NAs the vector size is equal to 0)
  if(sz2 == 0){
    return NA_REAL;
  }
  
  return *std::max_element(x2.begin(), x2.end()); //calculate the max value of the vector
}


//' Function that computes a maximum value for each bin
//'
//' @param x NumericVector
//' @param n intiger - number of bins
//' @keywords internal
// [[Rcpp::export]]
NumericVector binMax(NumericVector x,int n) {
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  
  // if the bins number larger than vector size, return zeros 
  if(sz < n){
    if(all(is_na(x))){  // if the vector consists of NAs, return NAs
      res=rep(NA_REAL,n);
    }
    return res;
  }
  
  double w_size = double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size,set the window size to 1
  if(sz == n){
    w_size = 1;
  }
  
  double prev = 0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2 ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    if(i == (n-1)){ // for the last bin
      end2 = sz;
    }
    prev = prev + w_size; // update the begining index of the slice
    
    NumericVector vec(&x[prev2], &x[end2]);
    res[i] = Max_c(vec);
  }
  
  return res;
}


//#-------------------------------------------------------------------------#
//' Function that computes a min value
//'
//' @param x NumericVector
//' @keywords internal
// [[Rcpp::export]]
double Min_c(NumericVector x){
  //remove NA from the vector
  NumericVector x2 = na_omit(x);
  
  int sz2 = x2.size() ;// get the length of the input vector after removing NAs
  
  //return NA when the vector consists of NAs (after removing NAs the vector size is equal to 0)
  if(sz2 == 0){
    return NA_REAL;
  }
  
  return *std::min_element(x2.begin(), x2.end()); //calculate the min value of the vector
}


//' Function that computes a minimum value for each bin
//'
//' @param x NumericVector
//' @param n intiger - number of bins
//' @keywords internal
// [[Rcpp::export]]
NumericVector binMin(NumericVector x,int n) {
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  
  // if the bins number larger than vector size, return zeros 
  if(sz < n){
    if(all(is_na(x))){  // if the vector consists of NAs, return NAs
      res=rep(NA_REAL,n);
    }
    return res;
  }
  
  double w_size = double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size,set the window size to 1
  if(sz == n){
    w_size = 1;
  }
  
  double prev = 0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    if(i == (n-1)){ // for the last bin
      end2 = sz;
    }
    prev = prev + w_size; // update the begining index of the slice
    NumericVector vec(&x[prev2], &x[end2]);
    res[i] = Min_c(vec);
  }
  
  return res;
}


//#-------------------------------------------------------------------------#
//' Function that computes a sum value
//'
//' @param x NumericVector
//' @keywords internal
// [[Rcpp::export]]
double Sum_c(NumericVector x){
  //remove NA from the vector
  NumericVector x2 = na_omit(x);
  
  int sz2 = x2.size() ;// get the length of the input vector after removing NAs
  
  //return NA when the vector consists of NAs (after removing NAs the vector size is equal to 0)
  if(sz2 == 0){
    return NA_REAL;
  }
  return std::accumulate(x2.begin(), x2.end(), 0.0); //calculate the sum value of the vector
}

//' Function that computes a sum of values in a bin
//'
//' @param x NumericVector - vector of values of a bin
//' @param n intiger - number of bins
//' @keywords internal
// [[Rcpp::export]]
NumericVector binSum(NumericVector x,int n) {
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  
  // if the bins number larger than vector size, return zeros 
  if(sz < n){
    if(all(is_na(x))){  // if the vector consists of NAs, return NAs
      res=rep(NA_REAL,n);
    }
    return res;
  }
  
  double w_size = double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size,set the window size to 1
  if(sz == n){
    w_size = 1;
  }
  
  double prev = 0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    if(i == (n-1)){ // for the last bin
      end2 = sz;
    }
    prev = prev + w_size; // update the begining index of the slice
    
    NumericVector vec(&x[prev2], &x[end2]);
    res[i] = Sum_c(vec);
  }
  
  return res;
}

//#-------------------------------------------------------------------------#
//' Function creates a matrix storing data with desirable number of bins for each window 
//'
//' listSliceMean() function calls the binMean() function
//'  
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix  listSliceMean(List xlist,int n) {
  int m = xlist.size();
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec = binMean(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}

//' Function creates a matrix storing data with desirable number of bins for each window 
//'
//' listSliceMean() function calls the binMedian() function 
//' 
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix  listSliceMedian(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec = binMedian(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}


//' Function creates a matrix storing data with desirable number of bins for each window 
//'
//' listSliceMax() function calls the binMax() function 
//' 
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix  listSliceMax(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec = binMax(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}


//' Function creates a matrix storing data with desirable number of bins for each window 
//'
//' listSliceMin() function calls the binMin() function 
//' 
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix  listSliceMin(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec = binMin(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}


//' Function creates a matrix storing data with desirable number of bins for each window 
//'
//' listSliceSum() function calls the binSum() function 
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix  listSliceSum(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec = binSum(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}

