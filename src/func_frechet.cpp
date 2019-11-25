/*
 * Though intrinsic mean is covered, we do not have the function to cover variation for that case. 
 * Therefore, in this script, all relevant computations will be covered for both ext/int 
 * mean and variation.
 */

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"
#include "RiemBaseExt_auxiliary.h"

//------------ AUXILIARY FUNCTIONS ------------------------------------------------
// (1) Distance
//    1-1. cppdist_1toN : distance from 1 observation to N samples

/*
* 1-1. cppdist_1toN : distance from 1 observation to N samples
*    mat  y         : matrix for 1 sample
*    cube X         : cube of N samples
*    string mfdname : manifold name
*    bool ext       : TRUE for extrinsic distance, FALSE for intrinsic
*/
arma::vec cppdist_1toN(arma::mat y, arma::cube X, std::string mfdname, bool ext){
  const int N = X.n_slices; // N slices
  arma::vec output(N,fill::zeros);
  arma::mat x;
  for (int i=0;i<N;i++){
    x = X.slice(i);
    if (ext){ // extrinsic distance
      output(i) = riemfunc_extdist(x, y, mfdname);
    } else {  // intrinsic distance
      output(i) = riemfunc_dist(x, y, mfdname);
    }
  }
  return(output);
}

/*
* 2-1. convert_3toM : from stacked data into a matrix of extrinsic column vectors
*    cube X         : stacked data cube
*    string mfdname : manifold name
* 
*/

//------------ MAIN FUNCTIONS ----------------------------------------------
// [[Rcpp::export]]
Rcpp::List frechet_extrinsic(arma::cube X, std::string mfdname){
  // preprocessing
  int N          = X.n_slices; // number of data
  arma::mat Xext = convert_3toM(X, mfdname); // columns are now samples 
  
  arma::mat testdata = X.slice(0); // size of the data
  int dnrow = testdata.n_rows; 
  int dncol = testdata.n_cols; 
  
  // compute the mean and projected mean 
  arma::vec extmean = arma::mean(Xext, 1); // find mean for each row.
  arma::mat inveqmu = riemfunc_invequiv(extmean,dnrow,dncol,mfdname); // projected mean
  arma::vec equivmu = riemfunc_equiv(inveqmu,dnrow,dncol,mfdname); 
  
  // extrinsic variation : (1/N)sum|extmean-data|^2 + |extmean-equivmu|^2
  double variation = 0.0;
  double tmp = 0.0;
  // ext1. mean squared error
  for (int i=0;i<N;i++){
    tmp = arma::as_scalar(arma::norm(extmean-Xext.col(i),"fro"));
    variation += (tmp*tmp);
  }
  variation /= static_cast<double>(N);
  // ext2. latter difference
  tmp = arma::as_scalar(arma::norm(extmean-equivmu,"fro"));
  variation += (tmp*tmp);
  
  Rcpp::List output;
  output["x"] = inveqmu;
  output["variation"] = variation;
  return(output);
}

// now we have intrinsic mean 'x' and all data 'X' in a cube
// [[Rcpp::export]]
arma::vec frechet_intrinsic(arma::mat x, arma::cube X, std::string mfdname){
  const int N = X.n_slices;
  arma::vec distvec(N,fill::zeros);
  arma::mat tmpdata;
  double tmpdist;
  for (int i=0;i<N;i++){
    tmpdata = X.slice(i);
    tmpdist = riemfunc_dist(x, tmpdata, mfdname);
    distvec(i) = tmpdist;
  }
  return(distvec);
}
