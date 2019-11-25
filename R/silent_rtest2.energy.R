#' Two-sample E-Test of Equal Distributions using Energy Statistics
#' 
#' @keywords internal
#' @noRd
rtest2.energy <- function(input1, input2, mc.iter=1234){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input1))!="riemdata"){
    stop("* rtest2.energy : the 'input1' must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  if ((class(input2))!="riemdata"){
    stop("* rtest2.energy : the 'input2' must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  # acquire manifold name
  mfdname = tolower(input1$name)
  # get data size
  m = length(input1$data)  # for class 1
  n = length(input2$data)  # for class 2

  #-------------------------------------------------------
  ## Prepare for using 'ENERGY' package
  # stack two datasets into one datacube
  dcube  = rclust_stack3d2(input1, input2)
  # compute distance
  dcdist = rclust_pdist_cube(dcube, mfdname, type="extrinsic")
  
  #-------------------------------------------------------
  ## Run the test
  output  = eqdist.etest(dcdist, distance=TRUE, sizes=c(m,n), R=mc.iter)
  return(output)
}

# Example with rvmf
# library(RiemBase)
# library(Directional)
# x = riemfactory(t(rvmf(100, c(1,0,0), 5)), name="sphere")
# y = riemfactory(t(rvmf(50, c(1,0,0), 5)), name="sphere")
# dxy = rtest2.energy(x, y, mc.iter = 100)
# dxy
