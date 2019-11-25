/*
 * (01) pdist_intrinsic : given cube, compute pairwise distance matrix in an intrinsic manner
 * (02) pdist_extrinsic : given cube, compute pairwise distance matrix in an extrinsic manner
 * (03) convert_3toM    : given cube, apply equivariant embedding and stack data as columns
 * 
 */



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"
#include "RiemBaseExt_auxiliary.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// (01) pdist_intrinsic : given cube, compute pairwise distance matrix in an intrinsic manner
// [[Rcpp::export]]
arma::mat pdist_intrinsic(arma::cube X, std::string mfdname){
  const int N = X.n_slices;
  arma::mat distmat(N,N,fill::zeros);
  
  arma::mat x;
  arma::mat y;
  double distval;
  for (int i=0;i<(N-1);i++){
    x = X.slice(i);
    for (int j=(i+1);j<N;j++){
      y = X.slice(j);
      distval = riemfunc_dist(x, y, mfdname);
      
      distmat(i,j) = distval;
      distmat(j,i) = distval;
    }
  }
  return(distmat);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// (02) pdist_extrinsic : given cube, compute pairwise distance matrix in an extrinsic manner
// [[Rcpp::export]]
arma::mat pdist_extrinsic(arma::cube X, std::string mfdname){
  const int N = X.n_slices;
  arma::mat distmat(N,N,fill::zeros);
  
  arma::mat x;
  arma::mat y;
  double distval;
  for (int i=0;i<(N-1);i++){
    x = X.slice(i);
    for (int j=(i+1);j<N;j++){
      y = X.slice(j);
      distval = riemfunc_extdist(x, y, mfdname);
      
      distmat(i,j) = distval;
      distmat(j,i) = distval;
    }
  }
  return(distmat);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// (03) convert_3toM    : given cube, apply equivariant embedding and stack data as columns
// [[Rcpp::export]]
arma::mat convert_3toM(arma::cube X, std::string mfdname){
  const int N = X.n_slices;
  // trial to capture length of a vector under equivariant embedding
  arma::mat testdata = X.slice(0);
  int dnrow = testdata.n_rows;
  int dncol = testdata.n_cols;
  arma::vec tester   = riemfunc_equiv(testdata,dnrow,dncol,mfdname);
  const int L = tester.n_elem; // length of an embedded vector
  
  // now convert into the matrix
  arma::mat output(L,N,fill::zeros);
  for (int i=0;i<N;i++){
    testdata = X.slice(i);
    output.col(i) = riemfunc_equiv(testdata,dnrow,dncol,mfdname);
  }
  return(output);
}
