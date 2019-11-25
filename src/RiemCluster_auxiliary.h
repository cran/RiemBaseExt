#ifndef RiemCluster_AUXILIARY_H
#define RiemCluster_AUXILIARY_H

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"

arma::mat cpp_equivariant(arma::cube X, std::string mfdname);
arma::mat cpp_distint(arma::cube X, std::string mfdname);
arma::mat cpp_distext(arma::cube X, std::string mfdname);
arma::mat cpp_dist2int(arma::cube X, arma::cube Y, std::string mfdname);
arma::mat cpp_dist2ext(arma::cube X, arma::cube Y, std::string mfdname);
arma::mat cpp_dispersion(arma::cube X, arma::mat xmean, std::string mfdname);

#endif
