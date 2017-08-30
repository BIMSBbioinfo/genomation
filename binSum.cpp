#include <Rcpp.h>
#include <math.h>

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
  NumericVector prev2(n) ;// integers for indices
  NumericVector end2(n)  ;
  double end;
  for(int i = 0; i <= n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2[i] = ceil(prev); // get the integer index for slices over vector
    end2[i] = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
  }
  // std::cout << prev2  << "\n" <<end2 << "\n";
  double partSum;
  NumericVector::iterator it;
  for(int i = 0; i < n; i++) {
    partSum = 0;
    if((n-1)-i){
   //caltulate the sum of values from the intervals excluding the last one
      for(it = &x[prev2[i]]; it != &x[end2[i]]; ++it) {
         partSum += *it;
      } 
     }else{  //calculate the sum of values from the last interval 
    for(it = &x[prev2[i]]; it <= &x[end2[i]]; ++it) {
         partSum += *it;
       }
  }
   res[i] = partSum/(&x[end2[i]]-&x[prev2[i]]); //calculate the mean value of the bin
  }
  
  return res;
}


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
  NumericVector prev2(n) ;// integers for indices
  NumericVector end2(n)  ;
  double end;
  for(int i = 0; i <= n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2[i] = ceil(prev); // get the integer index for slices over vector
    end2[i] = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
  }
  //std::cout << prev2  << "\n" <<end2 << "\n";
  for(int i = 0; i < n; i++) {
    int dint = &x[end2[i]] - &x[prev2[i]];  //size of the ith bin 
           //caltulate the median of values from the bin
      if(dint%2 == 0){
       std::sort(&x[prev2[i]], &x[end2[i]]);
     // double d = x[prev2[i]]+(dint/2)-1;  // (dint/2)
    //  double d2 = x[prev2[i]]+dint/2;  //  (dint/2)+1
      res[i] = (x[(prev2[i]+(dint/2)-1)] + (x[(prev2[i]+dint/2)]))/2;
     }else{
     std::sort(&x[prev2[i]], &x[end2[i]]);
     res[i] = x[(prev2[i] + dint/2)]; //(dint+1)/2
     }
    }
  return res;
}


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
  NumericVector prev2(n) ;// integers for indices
  NumericVector end2(n)  ;
  double end;
  for(int i = 0; i <= n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2[i] = ceil(prev); // get the integer index for slices over vector
    end2[i] = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
  }

  for(int i = 0; i < n; i++) {
    res[i] = *std::max_element(&x[prev2[i]], &x[end2[i]]); //calculate the max value in the bin
  }
  
  return res;
}


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
  NumericVector prev2(n) ;// integers for indices
  NumericVector end2(n)  ;
  double end;
  for(int i = 0; i <= n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2[i] = ceil(prev); // get the integer index for slices over vector
    end2[i] = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
  }
  
  for(int i = 0; i < n; i++) {
    res[i] = *std::min_element(&x[prev2[i]], &x[end2[i]]); //calculate the max value in the bin
  }
  return res;
}

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
  NumericVector prev2(n) ;// integers for indices
  NumericVector end2(n)  ;
  double end;
  for(int i = 0; i <= n; i++) {
    end = prev + (w_size); //get the end index of the interval
    prev2[i] = ceil(prev); // get the integer index for slices over vector
    end2[i] = ceil(end);
    prev = prev + w_size; // update the begining index of the slice
  }
  
  for(int i = 0; i < n; i++) {
    res[i] = std::accumulate(&x[prev2[i]], &x[end2[i]], 0.0); //calculate the max value in the bin
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix  listSliceMean(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  NumericVector tabx;
  for (int i = 0; i < m; i++) {
    subVec = binMean(xlist[i], n); //gives vector of mean values
    res(i, _) = subVec;             //adds the vector to the matrix
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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#x=binSum(1:35,10)
#x
#binMax(1:35,10)

#listSliceSum(list(1:35,1:25),10)
listSliceMean(list(1:35,1:25),10)
#listSliceMean3(list(1:35,1:25),10)
listSliceMean(list(c(120,1:35,120),1:25),10)
listSliceMedian(list(c(120,1:35,120),1:25),10)
listSliceMax(list(c(120,1:35,120),1:25),10)
listSliceMin(list(c(120,1:35,120),1:25),10)
listSliceSum(list(c(120,1:35,120),1:25),10)
*/
