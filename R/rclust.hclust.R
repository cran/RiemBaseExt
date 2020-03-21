#' Hierarchical Agglomerative Clustering for Manifold-valued Data
#' 
#' Hierarchical clustering is a generally applicable clustering algorithm 
#' as long as we have concept of dissimilarity. We adopt \code{hclust} algorithm 
#' by \pkg{fastcluster} package. See \code{\link[fastcluster]{hclust}} for more details.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link[RiemBase]{riemfactory}} for more details.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param method the agglomeration method to be used. This must be (an unambiguous abbreviation of) one of \code{"single"},
#' \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details.
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
#' ## compare extrinsic and intrinsic hierarchical clustering
#' hext <- rclust.hclust(data, type="extrinsic")
#' hint <- rclust.hclust(data, type="intrinsic")
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(hext, main="extrinsic+single")
#' plot(hint, main="intrinsic+single")
#' par(opar)
#' 
#' @references 
#' \insertRef{mullner_fastcluster_2013}{RiemBaseExt}
#' 
#' @export
rclust.hclust <- function(input, type=c("extrinsic","intrinsic"),
                          method=c("single","complete","average","mcquitty",
                                   "ward.D","ward.D2","centroid","median"),
                          members=NULL){
  #-------------------------------------------------------
  if (inherits(input, "dist")){
    dmethod = match.arg(method)
    output  = fastcluster::hclust(input, method=dmethod, members=members)
    return(output)
  } else {
    if (is.matrix(input)){
      input = RiemBase::riemfactory(t(input), name="euclidean")
    }
    # must be of 'riemdata' class
    if ((class(input))!="riemdata"){
      stop("* rclust.hclust : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
    }
    # acquire manifold name
    mfdname = tolower(input$name)
    # type and method checking checking
    mfdtype = match.arg(type)
    dmethod = match.arg(method)
    
    #-------------------------------------------------------
    # compute distance and apply 'fastcluster::hclust'
    dist   = stats::as.dist(rclust_pdist(input, type=mfdtype))
    output = fastcluster::hclust(dist, method=dmethod, members=members)
    return(output) 
  }
}



# # Old Example with 'sphere3d' data ---------------------------------------- 
# ## load point clouds of 3 groups on sphere
# data(sphere3)
# ## visualize data on sphere
# library("rgl")
# mfrow3d(1,2)
# spheres3d(0,0,0,lit=FALSE,color="white")
# spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")
# spheres3d(sphere3[,1],sphere3[,2],sphere3[,3],col="black",radius=0.03)
# 
# ## wrap the data into 'riemdata' class
# ## note that when given as a matrix, each data be a column.
# library(RiemBase)
# rsp3 <- riemfactory(t(sphere3), name="Sphere")
# 
# ## apply hierarchical clustering with complete linkage
# result  <- rclust.hclust(rsp3, method="complete")
# 
# ## get clustering labels under 3-group assumption
# label3  <- stats::cutree(result, k=3)
# 
# ## find indices corresponding to 3 groups
# idx.gp1 <- which(label3==unique(label3)[1])
# idx.gp2 <- which(label3==unique(label3)[2])
# idx.gp3 <- which(label3==unique(label3)[3])
# 
# ## visualize clustering results of 3 groups
# rgl::next3d()
# spheres3d(0,0,0,lit=FALSE,color="white")
# spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")
# spheres3d(sphere3[idx.gp1,1],sphere3[idx.gp1,2],sphere3[idx.gp1,3],col="red",radius=0.03)
# spheres3d(sphere3[idx.gp2,1],sphere3[idx.gp2,2],sphere3[idx.gp2,3],col="blue",radius=0.03)
# spheres3d(sphere3[idx.gp3,1],sphere3[idx.gp3,2],sphere3[idx.gp3,3],col="yellow",radius=0.03)
# 
# ## you can close the RGL figure by 'rgl.close()'

