#' Pairwise Distance for Two Sets of Data on Manifolds
#' 
#' For two points \eqn{x,y \in \mathcal{M}}, two modes of distances are available; \emph{intrinsic} for geodesic distance 
#' on the manifold and \emph{extrinsic} for standard norm after equivariant embedding into Euclidean space. This function 
#' differs from \code{\link{rstat.pdist}} in that it now computes distances between two sets of data.
#' 
#' @param input1 a S3 object of \code{riemdata} class of \eqn{N} objects.
#' @param input2 a S3 object of \code{riemdata} class of \eqn{M} objects. See \code{\link[RiemBase]{riemfactory}} for more details.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' 
#' @return an \eqn{(N\times M)} matrix of pairwise distances.
#' 
#' @seealso \code{\link{rstat.pdist}}
#' 
#' @examples
#' ### Generate 100 data points on Sphere S^2.
#' #   50 points from near (0,0,1), and
#' #   50 points from near (0,0,-1)
#' 
#' ndata = 50
#' theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
#' tmpx  = cos(theta)
#' tmpy  = sin(theta)
#' 
#' ### Wrap those as 'riemdata' class
#' data1 = list()
#' data2 = list()
#' for (i in 1:ndata){
#'   tgt1 = c(tmpx[i],tmpy[i],1)  + stats::rnorm(3,sd=0.1)
#'   tgt2 = c(tmpx[i],tmpy[i],-1) + stats::rnorm(3,sd=0.1)
#'   
#'   data1[[i]] = tgt1/sqrt(sum(tgt1^2)) # projection near (0,0,1)
#'   data2[[i]] = tgt2/sqrt(sum(tgt2^2)) #                 (0,0,-1)  
#' }
#' spdata1 = RiemBase::riemfactory(data1, name="sphere")
#' spdata2 = RiemBase::riemfactory(data2, name="sphere")
#' 
#' ### Compute Two Types of Distances and Visualize
#' dist.int = rstat.pdist2(spdata1, spdata2, type="intrinsic")
#' dist.ext = rstat.pdist2(spdata1, spdata2, type="extrinsic")
#' 
#' ### Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' image(dist.int, main="intrinsic")
#' image(dist.ext, main="extrinsic")
#' par(opar)
#' 
#' @export
rstat.pdist2 <- function(input1, input2, type=c("intrinsic","extrinsic")){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if (!inherits(input1,"riemdata")){
    stop("* rstat.pdist2 : the 'input1' must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  if (!inherits(input2,"riemdata")){
    stop("* rstat.pdist2 : the 'input2' must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  if (input1$name != input2$name){
    stop("* rstat.pdist2 : two inputs should be of same manifold.")
  }
  # type checking
  mytype  = match.arg(type)

  #-------------------------------------------------------
  # compute
  outdist = rclust_pdist2(input1, input2, type=mytype)
  return(outdist)
}
