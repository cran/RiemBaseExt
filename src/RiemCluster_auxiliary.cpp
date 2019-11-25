// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"
#include "RiemCluster_auxiliary.h"


/*
* (01) cpp_equivariant : convert manifold-valued data into embedded ones
* (02) cpp_distint     : intrinsic distance between pairs of datapoints
*      cpp_dist2int    : dist2 version
* (03) cpp_distext     : extrinsic distance
*      cpp_dist2ext    : dist2 version
* (04) cpp_dispersion  : given a mean and targets, compute log-pulled vectors in TpM as columns
*/

////////////////////////////////////////////////////////////////////////
// (01) cpp_equivariant : convert manifold-valued data into embedded ones.
// [[Rcpp::export]]
arma::mat cpp_equivariant(arma::cube X, std::string mfdname){
  // taking one exempler data
  arma::mat testdata = X.slice(0);
  // get some relevant integer parameters
  int N = X.n_slices;
  int dnrow = testdata.n_rows;
  int dncol = testdata.n_cols;
  // apply equivariant embedding for the testdata to get the size.
  arma::vec testequiv = riemfunc_equiv(testdata,dnrow,dncol,mfdname);
  int P = testequiv.n_elem;
  
  // now apply for each data, stacked as rows
  arma::mat output(N,P,fill::zeros);
  for (int i=0;i<N;i++){
    testdata  = X.slice(i);
    testequiv = riemfunc_equiv(testdata,dnrow,dncol,mfdname);
    output.row(i) = testequiv.t();
  }
  return(output);
}

////////////////////////////////////////////////////////////////////////
// (02) cpp_distint     : intrinsic distance between pairs of datapoints
// [[Rcpp::export]]
arma::mat cpp_distint(arma::cube X, std::string mfdname){
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
// [[Rcpp::export]]
arma::mat cpp_dist2int(arma::cube X, arma::cube Y, std::string mfdname){
  const int N = X.n_slices;
  const int M = Y.n_slices;
  arma::mat distmat(N,M,fill::zeros);
  
  arma::mat x;
  arma::mat y;
  for (int i=0;i<N;i++){
    x = X.slice(i);
    for (int j=0;j<M;j++){
      y = Y.slice(j);
      distmat(i,j) = riemfunc_dist(x,y,mfdname);
    }
  }
  return(distmat);
}
////////////////////////////////////////////////////////////////////////
// (03) cpp_distext     : extrinsic distance
// [[Rcpp::export]]
arma::mat cpp_distext(arma::cube X, std::string mfdname){
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
// [[Rcpp::export]]
arma::mat cpp_dist2ext(arma::cube X, arma::cube Y, std::string mfdname){
  const int N = X.n_slices;
  const int M = Y.n_slices;
  arma::mat distmat(N,M,fill::zeros);
  
  arma::mat x;
  arma::mat y;
  for (int i=0;i<N;i++){
    x = X.slice(i);
    for (int j=0;j<M;j++){
      y = Y.slice(j);
      distmat(i,j) = riemfunc_extdist(x,y,mfdname);
    }
  }
  return(distmat);
}

////////////////////////////////////////////////////////////////////////
// (04) cpp_dispersion  : given a mean and targets, compute log-pulled vectors in TpM as columns
// [[Rcpp::export]]
arma::mat cpp_dispersion(arma::cube X, arma::mat xmean, std::string mfdname){
  // get some parameters
  int dnrow = xmean.n_rows;
  int dncol = xmean.n_cols;
  int nsubj = X.n_slices;
  
  // record vectorized logmaps
  arma::mat logvecs(dnrow*dncol, nsubj, fill::zeros);
  for (int i=0;i<nsubj;i++){
    logvecs.col(i) = arma::vectorise(riemfunc_log(xmean, X.slice(i), mfdname));
  }
  return(logvecs);
}
