#' DBSCAN for Manifold-valued Data
#' 
#' Density-Based Spatial Clustering of Applications with Noise (DBSCAN) is a generally applicable clustering algorithm 
#' as long as we have concept of dissimilarity. We adopt \code{dbscan} algorithm 
#' by \pkg{dbscan} package. See \code{\link[dbscan]{dbscan}} for more details.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link[RiemBase]{riemfactory}} for more details.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param eps size of the epsilon neighborhood.
#' @param minPts number of minimum points in the eps region (for core points). Default is 5 points.
#' @param ... extra parameters for DBSCAN algorithms including \code{eps} or \code{minPts}. See \code{\link[dbscan]{dbscan}} for more details.
#' 
#' @return an object of class \code{dbscan_fast} containing
#' \describe{
#'   \item{eps}{value of the eps parameter.}
#'   \item{minPts}{value of the minPts parameter.}
#'   \item{cluster}{An integer vector with cluster assignments. Zero indicates noise points.}
#' }
#' 
#' @examples
#' ## generate 50 points near (0,0,1)  and
#' ##          50 points near (0,0,-1) on Sphere S^2 
#' ndata = 50
#' theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
#' tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
#' tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
#' 
#' ## wrap it as 'riemdata' class
#' data  = list()
#' for (i in 1:ndata){
#'   tgt = c(tmpx[i],tmpy[i],1)
#'   data[[i]] = tgt/sqrt(sum(tgt^2)) # project onto Sphere
#' }
#' for (i in 1:ndata){
#'   tgt = c(tmpx[i],tmpy[i],-1)
#'   data[[i+ndata]] = tgt/sqrt(sum(tgt^2)) # project onto Sphere
#' }
#' data = RiemBase::riemfactory(data, name="sphere")
#' 
#' ## compare extrinsic and intrinsic DBSCAN
#' dbext <- rclust.dbscan(data, eps=0.5, type="extrinsic")
#' dbint <- rclust.dbscan(data, eps=0.5, type="intrinsic")
#' 
#' ## let's visualize the results via MDS in R^2
#' pdist = stats::as.dist(RiemBase::rbase.pdist(data))
#' dat2d = stats::cmdscale(pdist, k=2)
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(dat2d[,1], dat2d[,2], col=dbext$cluster, pch=19, main="extrinsic+dbscan")
#' plot(dat2d[,1], dat2d[,2], col=dbint$cluster, pch=19, main="intrinsic+dbscan")
#' par(opar)
#' 
#' @seealso See \code{\link[dbscan]{dbscan}} for details.
#' 
#' @references 
#' \insertRef{ester_density-based_1996}{RiemBaseExt}
#' 
#' \insertRef{hahsler_dbscan_2019}{RiemBaseExt}
#' 
#' @export
rclust.dbscan <- function(input, type=c("extrinsic","intrinsic"), eps, minPts=5, ...){
  #-------------------------------------------------------
  if (inherits(input, "dist")){  # METHOD 1 : dist object
    output  = dbscan::dbscan(input, eps, minPts=minPts, ...)
    return(output)
  } else {                       # METHOD 2 : we need to compute distance
    if (is.matrix(input)){
      input = RiemBase::riemfactory(t(input), name="euclidean")
    }
    # must be of 'riemdata' class
    if ((class(input))!="riemdata"){
      stop("* rclust.dbscan : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
    }
    # acquire manifold name
    mfdname = tolower(input$name)
    # type and method checking checking
    mfdtype = match.arg(type)

    #-------------------------------------------------------
    # compute distance and apply 'dbscan::dbscan'
    distobj = stats::as.dist(rclust_pdist(input, type=mfdtype))
    output  = dbscan::dbscan(distobj, eps, minPts=minPts, ...)
    return(output)
  }
}
