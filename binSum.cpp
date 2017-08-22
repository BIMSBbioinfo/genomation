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

// function to test max min and indeces shouldn't be exported
// [[Rcpp::export]]
double basic_function(NumericVector x,int prev,int end) { // this function is just a trial, can be deleted
  double res = *std::max_element(x.begin()+prev,x.begin()+end);// get max
  return res;
}


// [[Rcpp::export]]
std::vector<double>  binMax(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  
  std::vector<double>  res(n);// create the output vector
  
  double prev=0; // index for start positions over vector
  double w_size=double(sz-1)/double(n); // window size can be a double
  //std::cout << w_size << " || \n";
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz==n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  //w_size=round(w_size);
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n-1; i++) {
    end=prev+(w_size); //get the end index of vector
    
    prev2= ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    std::cout << prev2  << " " <<end2 << "\n";
    
    res[i] = *std::max_element(x.begin()+prev2,x.begin()+end2);// get max
    prev=prev+ w_size; // update the begining index of the slice
  }
  
  // for the last slice do calc outside the loop to be able to include last element in the slice
  // there might be more elegant way to do this but there is no time
  end=prev+(w_size); 
  //std::cout << prev  << " " <<end << "\n";
  prev2= ceil(prev);
  end2 = ceil(end);
  //std::cout << prev2  << " " <<end2 << "\n";
  res[(n-1)] = *std::max_element(x.begin()+prev2,x.end());
  //std::cout <<  *std::max(x.begin()+prev2,x.end()) << "\n";
  return res;
}

// function that calculates the bin sum
std::vector<double>  binSum(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  
  std::vector<double>  res(n);// create the output vector
  
  double prev=0; // index for start positions over vector
  double w_size=double(sz-1)/double(n); // window size can be a double
  //std::cout << w_size << " || \n";
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz==n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  //w_size=round(w_size);
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n-1; i++) {
    end=prev+(w_size); //get the end index of vector
    //std::cout << prev  << " " <<end << "\n";
    
    prev2= ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    res[i] = std::accumulate(x.begin()+prev2,x.begin()+end2, 0.0);// sum up
    prev=prev+ w_size; // update the begining index of the slice
  }
  
  // for the last slice do calc outside the loop to be able to include last element in the slice
  // there might be more elegant way to do this but there is no time
  end=prev+(w_size); 
  //std::cout << prev  << " " <<end << "\n";
  prev2= ceil(prev);
  end2 = ceil(end);
  //std::cout << prev2  << " " <<end2 << "\n";
  res[n-1] = std::accumulate(x.begin()+prev2,x.begin()+(end2+1), 0.0);
  return res;
}

std::vector<double>  binMean(NumericVector x,int n) {
  
  int sz = x.size() ;// get the length of the input vector
  
  std::vector<double>  res(n);// create the output vector
  
  double prev=0; // index for start positions over vector
  double w_size=double(sz-1)/double(n); // window size can be a double
  //std::cout << w_size << " || \n";
  
  // if the bins equals the vector size ,set the window size to 1
  if(sz==n){
    w_size=1;
  }
  
  // if the bins number larger than vector size return zeros 
  if(sz < n){
    return res;
  }
  
  //w_size=round(w_size);
  int prev2 ;// integers for indices
  int end2  ;
  double end;
  for(int i = 0; i < n-1; i++) {
    end=prev+(w_size); //get the end index of vector
    //std::cout << prev  << " " <<end << "\n";
    
    prev2= ceil(prev); // get the integer index for slices over vector
    end2 = ceil(end);
    res[i] = std::accumulate(x.begin()+prev2,x.begin()+end2, 0.0);// sum up
    res[i] = res[i]/(end2-prev2); // get mean
    prev=prev+ w_size; // update the begining index of the slice
  }
  
  // for the last slice do calc outside the loop to be able to include las slice
  // there might be more elegant way to do this but there is no time
  end=prev+(w_size); 
  //std::cout << prev  << " " <<end << "\n";
  prev2= ceil(prev);
  end2 = ceil(end);
  //std::cout << prev2  << " " <<end2 << "\n";
  res[n-1] = std::accumulate(x.begin()+prev2,x.begin()+(end2+1), 0.0)/(end2-prev2+1);
  return res;
}





// [[Rcpp::export]]
NumericMatrix  listSliceSum(List xlist,int n) {
  int m = xlist.size(); 
  std::cout << m  << " " <<n << "\n";
  NumericMatrix res(m, n);
  std::vector<double>  subVec;
  for (int i = 0; i < m; i++) {
    subVec=binSum(xlist[i],n);
    for (int j = 0; j < n; j++) {
      res(i, j)=subVec[j];
    }
  }
  
  return res;
}


// [[Rcpp::export]]
NumericMatrix  listSliceMean(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  std::vector<double>  subVec;
  for (int i = 0; i < m; i++) {
    subVec=binMean(xlist[i],n);
    for (int j = 0; j < n; j++) {
      res(i, j)=subVec[j];
    }
  }
  
  return res;
}

// [[Rcpp::export]]
NumericVector binMean2(NumericVector x,int n) {
  int p = x.size();  // get the length of the input vector
  double w_size=double(p-1)/double(n); // window size can be a double
  int step = ceil(w_size);
  NumericVector res(n);// create the output vector
  
  NumericMatrix bin(step,n);
  
  for (int i = 0; i < p; i++){
    bin[i] = x[i];
  }
  
  NumericVector bz(step);
  
  for(int z = 0; z < n; z++){
    bz = bin( _, z);
    res[z] = std::accumulate(bz.begin(), bz.end(), 0.0)/bz.size();
  } 
  
  return res;
}



// [[Rcpp::export]]
NumericMatrix  listSliceMean2(List xlist,int n) {
  int m = xlist.size(); 
  NumericMatrix res(m, n);
  NumericVector  subVec;
  for (int i = 0; i < m; i++) {
    subVec = binMean2(xlist[i],n);
    res(i, _) = subVec;
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
binMax(1:35,10)

listSliceSum(list(1:35,1:25),10)
listSliceMean(list(1:35,1:25),10)
listSliceMean2(list(1:35,1:25),10)
listSliceMean(list(c(120,1:35,120),1:25),10)
listSliceMean2(list(c(120,1:35,120),1:25),10)
*/
