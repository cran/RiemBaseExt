#' Fréchet Mean and Variation of Manifold-valued Data
#' 
#' For manifold-valued data, Fréchet mean is the solution of following cost function,
#' \deqn{\textrm{min}_x \sum_{i=1}^n \rho^2 (x, x_i),\quad x\in\mathcal{M}}
#' for a given data \eqn{\{x_i\}_{i=1}^n} and \eqn{\rho(x,y)} is the geodesic distance 
#' between two points on manifold \eqn{\mathcal{M}}. It uses a gradient descent method 
#' with a backtracking search rule for updating.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link[RiemBase]{riemfactory}} for more details.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param int.eps stopping criterion for the norm of gradient.
#' @param parallel a flag for enabling parallel computation.
#' 
#' @return a named list containing
#' \describe{
#' \item{mu}{an estimated Fréchet mean matrix.}
#' \item{variation}{Fréchet variation with the estimated mean.}
#' }
#' 
#' @examples
#' ### Generate 50 data points on Sphere S^2 near (0,0,1).
#' ndata = 50
#' theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
#' tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
#' tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
#' 
#' ### Wrap it as 'riemdata' class
#' data  = list()
#' for (i in 1:ndata){
#'   tgt = c(tmpx[i],tmpy[i],1)
#'   data[[i]] = tgt/sqrt(sum(tgt^2)) # project onto Sphere
#' }
#' data = RiemBase::riemfactory(data, name="sphere")
#' 
#' ### Compute Fréchet Mean and Variation
#' out1 = rstat.frechet(data)                   # intrinsic
#' out2 = rstat.frechet(data,parallel=TRUE)     # parallel implementation
#' out3 = rstat.frechet(data, type="extrinsic") # extrinsic
#' 
#' @references 
#' \insertRef{karcher_riemannian_1977}{RiemBaseExt}
#' 
#' \insertRef{kendall_probability_1990}{RiemBaseExt}
#' 
#' \insertRef{afsari_convergence_2013}{RiemBaseExt}
#' 
#' @export
rstat.frechet <- function(input, type=c("intrinsic","extrinsic"), int.eps=1e-6, parallel=FALSE){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input))!="riemdata"){
    stop("* rstat.frechet : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  # acquire manifold name
  mfdname = tolower(input$name)
  # type checking
  type = match.arg(type)
  # parallel flag
  pflag = parallel
  if (!is.logical(pflag)){
    stop("* rstat.frechet : 'parallel' must be a logical flag.")
  }
  
  #-------------------------------------------------------
  # TYPE BRANCHING : Extrinsic and Intrinsic
  newdata = aux_stack3d(input) # 3d-stacked data
  if (type=="intrinsic"){
    tmpout  = RiemBase::rbase.mean(input, eps=int.eps, parallel = pflag)
    tmpsol  = tmpout$x
    distvec = frechet_intrinsic(tmpsol, newdata, mfdname);
    
    output = list()
    output$mu = tmpsol
    output$variation = (sum(distvec^2))/(dim(newdata)[3])
  } else {
    result = frechet_extrinsic(newdata, mfdname)
    output = list()
    output$mu        = result$x
    output$variation = result$variation
  }
  
  
  #-------------------------------------------------------
  # RETURN
  return(output)
}


#' @keywords internal
#' @noRd
rstat.frechet.cube <- function(datacube, mfdname, type=c("extrinsic","intrinsic"), int.eps=1e-6, parallel=FALSE){
  #-------------------------------------------------------
  # type checking
  type = match.arg(type)
  # parallel flag
  pflag = parallel
  if (!is.logical(pflag)){
    stop("* rstat.frechet : 'parallel' must be a logical flag.")
  }
  
  #-------------------------------------------------------
  # TYPE BRANCHING : Extrinsic and Intrinsic
  newdata = datacube # 3d-stacked data
  if (type=="intrinsic"){
    mean.internal = utils::getFromNamespace("rbase.mean.cube","RiemBase")
    tmpout = mean.internal(datacube, mfdname, eps=int.eps, parallel=pflag)
    tmpsol  = tmpout$x
    distvec = frechet_intrinsic(tmpsol, newdata, mfdname);
    
    output = list()
    output$x = tmpsol
    output$variation = (sum(distvec^2))/(dim(newdata)[3])
  } else {
    output  = frechet_extrinsic(newdata, mfdname);
  }
  
  
  #-------------------------------------------------------
  # RETURN
  return(output)
}