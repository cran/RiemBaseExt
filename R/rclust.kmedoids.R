#' k-Medoids Clustering Manifold-valued Data
#' 
#' k-Medoids is a generally applicable clustering algorithm 
#' as long as we have concept of dissimilarity. We adopt \code{pam} algorithm 
#' by \pkg{cluster} package. See \code{\link[cluster]{pam}} for more details.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link{riemfactory}} for more details.
#' @param k the number of clusters.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' 
#' @return an object of class \code{pam}. See \code{\link[cluster]{pam}} for details.
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
#' ## compare k-medoids with intrinsic distance for different k's
#' k2 <- rclust.kmedoids(data, k=2, type="intrinsic")
#' k3 <- rclust.kmedoids(data, k=3, type="intrinsic")
#' k4 <- rclust.kmedoids(data, k=4, type="intrinsic")
#' 
#' ## let's visualize the results via MDS in R^2
#' pdist = stats::as.dist(RiemBase::rbase.pdist(data))
#' dat2d = stats::cmdscale(pdist, k=2)
#' 
#' opar <- par(mfrow=c(1,3))
#' plot(dat2d[,1], dat2d[,2], col=k2$clustering, pch=19, main="k=2")
#' plot(dat2d[,1], dat2d[,2], col=k3$clustering, pch=19, main="k=3")
#' plot(dat2d[,1], dat2d[,2], col=k4$clustering, pch=19, main="k=4")
#' par(opar)
#' 
#' @export
rclust.kmedoids <- function(input, k=2, type=c("extrinsic","intrinsic")){
  #-------------------------------------------------------
  myk = round(k)
  if (inherits(input, "dist")){
    output  = cluster::pam(input, k=myk)
    return(output)
  } else {
    if (is.matrix(input)){
      input = RiemBase::riemfactory(t(input), name="euclidean")
    }
    # must be of 'riemdata' class
    if ((class(input))!="riemdata"){
      stop("* rclust.kmedoids : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
    }
    # acquire manifold name
    mfdname = tolower(input$name)
    # type and method checking checking
    mfdtype = match.arg(type)

    #-------------------------------------------------------
    # compute distance and apply 'fastcluster::hclust'
    distobj = stats::as.dist(rclust_pdist(input, type=mfdtype))
    output  = cluster::pam(distobj, k=myk)
    return(output) 
  }
}