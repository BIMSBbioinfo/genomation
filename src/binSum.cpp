#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

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
    
    int dint = &x[end2] - &x[prev2];  //size of the ith bin 
    //caltulate the median of values from the bin
    if(dint%2 == 0){
      std::sort(&x[prev2], &x[end2]);
      // double d = x[prev2]+(dint/2)-1;  // (dint/2)
      //  double d2 = x[prev2]+dint/2;  //  (dint/2)+1
      res[i] = (x[(prev2+(dint/2)-1)] + (x[(prev2+dint/2)]))/2;
    }else{
      std::sort(&x[prev2], &x[end2]);
      res[i] = x[(prev2 + dint/2)]; //(dint+1)/2
    }
    
  }
  
  return res;
}

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

// [[Rcpp::export]]
NumericMatrix  listSliceMean(List xlist,int n) {
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

// [[Rcpp::export]]
NumericMatrix  listSliceMedian(List xlist,int n) {
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

// [[Rcpp::export]]
NumericMatrix  listSliceMax(List xlist,int n) {
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


// [[Rcpp::export]]
NumericMatrix  listSliceMin(List xlist,int n) {
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


// [[Rcpp::export]]
NumericMatrix  listSliceSum(List xlist,int n) {
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
