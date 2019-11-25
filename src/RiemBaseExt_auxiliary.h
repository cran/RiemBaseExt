#ifndef RiemBaseExt_AUXILIARY_H
#define RiemBaseExt_AUXILIARY_H

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"

arma::mat pdist_intrinsic(arma::cube X, std::string mfdname);
arma::mat pdist_extrinsic(arma::cube X, std::string mfdname);
arma::mat convert_3toM(arma::cube X, std::string mfdname);

#endif
