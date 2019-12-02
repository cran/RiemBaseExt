#' Extension to 'RiemBase' Package 
#' 
#' We provide a number of off-the-shelf algorithms for clustering, elementary computation, and others. 
#' See Bhattacharya and Bhattacharya (2012) <doi:10.1017/CBO9781139094764> for more details regarding Statistics on Manifolds.
#' 
#' @docType package
#' @name RiemBaseExt-package
#' @import Rdpack
#' @importFrom cluster pam
#' @importFrom dbscan dbscan
#' @importFrom energy eqdist.etest
#' @importFrom fastcluster hclust
#' @importFrom kernlab specc
#' @importFrom parallel detectCores
#' @importFrom stats as.dist cutree sd var na.omit
#' @importFrom utils getFromNamespace head tail packageVersion
#' @importFrom RiemBase rbase.pdist riemfactory
#' @importFrom Rcpp evalCpp
#' @useDynLib RiemBaseExt
NULL