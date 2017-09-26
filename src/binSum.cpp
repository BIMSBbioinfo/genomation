#include <math.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
//#-------------------------------------------------------------------------#
//' Functions that compute the value for each bin
//'
//' binMean() computes a mean value, binMedian() computes a median value, binMax() and binMin() 
//' give maximum and minumum values respectively, binSum() computes the sum of the values in a bin
//'
//' @param x NumericVector - vector of values of a bin
//' @param n intiger - number of bins
//' @export
//' @rdname bin
// [[Rcpp::export]]
NumericVector binMean(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz == n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  double prev=0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i< n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
    res[i] = std::accumulate(&x[prev2], &x[end2], 0.0)/(&x[end2]-&x[prev2]); //calculate the mean value of the bin
  }
  
  return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
double Median_c(NumericVector x){
  int dint = x.size();
  double res;
  if(dint%2 == 0){
    std::sort(x.begin(), x.end());
    res = (x[(dint/2)-1] + x[dint/2] ) / 2; 
  }else{
    std::nth_element(x.begin(), x.begin()+dint/2, x.end());
    res = x[dint/2];
  }
  return res;
}

//' @export
//' @rdname bin
// [[Rcpp::export]]
NumericVector binMedian(NumericVector x, int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz == n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  double prev=0; // index for start positions over vector
  int prev2 ;  // integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
    
    NumericVector vec(&x[prev2], &x[end2]);
    res[i] = Median_c(vec);
    
  }
  
  return res;
}


//' @export
//' @rdname bin
// [[Rcpp::export]]
NumericVector binMax(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz == n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  double prev=0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2 ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
    res[i] = *std::max_element(&x[prev2], &x[end2]); //calculate the max value in the bin
  }
  
  return res;
}

//' @export
//' @rdname bin
// [[Rcpp::export]]
NumericVector binMin(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz == n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  double prev=0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
    res[i] = *std::min_element(&x[prev2], &x[end2]); //calculate the max value in the bin
  }
  
  return res;
}

//' @export
//' @rdname bin
// [[Rcpp::export]]
NumericVector binSum(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  NumericVector res(n);// create the output vector
  double w_size=double(sz)/double(n); // window size can be a double
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz == n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  double prev=0; // index for start positions over vector
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2 = ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
    res[i] = std::accumulate(&x[prev2], &x[end2], 0.0); //calculate the sum value of the bin
    
  }
  
  return res;
}


//#-------------------------------------------------------------------------#
//' Function reverses a vector 
//' 
//' @param y NumericVector
//' @export
//' @rdname reverse_L
// [[Rcpp::export]]
NumericVector reverse_L(NumericVector y) {
  std::reverse(y.begin(), y.end());
  return y;
}

//#-------------------------------------------------------------------------#
//' Functions create a matrix storing the data with desirable number of bins for each window 
//'
//' listSliceMean() calls the binMean() function, listSliceMedian calls the binMedian(), 
//' listSliceMax() - binMax(), listSliceMin() - binMin(), listSliceSum() - binSum()
//'
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @param ranks CharacterVector - position of the windows whose strand is "-"
//' @export
//' @rdname listSlice
// [[Rcpp::export]]
NumericMatrix  listSliceMean(List xlist,int n, CharacterVector ranks)   {
  if(ranks.length()){
    CharacterVector f = xlist.names();
    std::string p;
    for(int j=0; j < ranks.size(); j++){
      p = *std::find(f.begin(),f.end(), ranks(j));
      xlist[p] = reverse_L(xlist[p]);
    }
  }
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec=binMean(xlist[i],n);
    for (int j = 0; j < n; j++) {
      res(i, j)=subVec[j];
    }
  }
  return res;
}

//' @export
//' @rdname listSlice
// [[Rcpp::export]]
NumericMatrix  listSliceMedian(List xlist,int n, CharacterVector ranks) {
  if(ranks.length()){
    CharacterVector f = xlist.names();
    std::string p;
    for(int j=0; j < ranks.size(); j++){
      p = *std::find(f.begin(),f.end(), ranks(j));
      xlist[p] = reverse_L(xlist[p]);
    }
  }
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  NumericVector tabx;
  for (int i = 0; i < m; i++) {
    subVec = binMedian(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}

//' @export
//' @rdname listSlice
// [[Rcpp::export]]
NumericMatrix  listSliceMax(List xlist,int n, CharacterVector ranks) {
  if(ranks.length()){
    CharacterVector f = xlist.names();
    std::string p;
    for(int j=0; j < ranks.size(); j++){
      p = *std::find(f.begin(),f.end(), ranks(j));
      xlist[p] = reverse_L(xlist[p]);
    }
  }
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  NumericVector tabx;
  for (int i = 0; i < m; i++) {
    subVec = binMax(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}

//' @export
//' @rdname listSlice
// [[Rcpp::export]]
NumericMatrix  listSliceMin(List xlist,int n, CharacterVector ranks) {
  if(ranks.length()){
    CharacterVector f = xlist.names();
    std::string p;
    for(int j=0; j < ranks.size(); j++){
      p = *std::find(f.begin(),f.end(), ranks(j));
      xlist[p] = reverse_L(xlist[p]);
    }
  }
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  NumericVector tabx;
  for (int i = 0; i < m; i++) {
    subVec = binMin(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}

//' @export
//' @rdname listSlice
// [[Rcpp::export]]
NumericMatrix  listSliceSum(List xlist,int n, CharacterVector ranks) {
  if(ranks.length()){
    CharacterVector f = xlist.names();
    std::string p;
    for(int j=0; j < ranks.size(); j++){
      p = *std::find(f.begin(),f.end(), ranks(j));
      xlist[p] = reverse_L(xlist[p]);
    }
  }
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  NumericVector tabx;
  for (int i = 0; i < m; i++) {
    subVec = binSum(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
  }
  return res;
}


//#-------------------------------------------------------------------------#
//' Function reorders matrix to obtain the original order of the windows
//'
//' @param x NumericMatrix 
//' @param p NumericVector - stors an original window order
//' @export 
//' @rdname ranksOrder
// [[Rcpp::export]]
NumericMatrix  ranksOrder(NumericMatrix x, NumericVector p) {
  int m = x.nrow();
  int n = x.ncol();
  
  NumericMatrix res(m, n);
  NumericVector r = clone(p);
  
  std::sort(r.begin(), r.end());
  
  for (int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      if(p[i]==r[j]){
        res.row(j) = x.row(i);
        break;
      }
    }
  }
  return res;
}


//#-------------------------------------------------------------------------#
//' Function computes a matrix that stores the data with desirable number of bins and 
//' keeps original order of the windows  
//'
//' @param xlist List of vectors storing values of a bin
//' @param n intiger - number of bins
//' @param negranks CharacterVector - position of the windows whose strand is "-"
//' @param binOp - "bin" option - "mean", "max", "sum", "median", "min"
//' @param p NumericVector - stors an original window order
//' @export 
//' @rdname matRes
// [[Rcpp::export]]
NumericMatrix matRes(List xlist, int n, CharacterVector negranks, std::string binOp, NumericVector p){
  NumericMatrix res(xlist.size(),n);
  
 if(binOp.compare("mean") == 0){
   res = listSliceMean(xlist, n, negranks);
 }else if(binOp.compare("max") == 0){
   res = listSliceMax(xlist, n, negranks);
 }else if(binOp.compare("sum") == 0){
   res = listSliceSum(xlist, n, negranks);
 }else if(binOp.compare("median") == 0){
   res = listSliceMedian(xlist, n, negranks);
 }else if(binOp.compare("min") == 0){
   res = listSliceMin(xlist, n, negranks); 
 }
 
  return ranksOrder(res, p);
    
}

