#' Spectral Clustering : only extrinsic is valid for kernel issue.
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' those need to be found and added to Zotero
#' 
#' @references 
#' \insertRef{kernlab}{RiemBaseExt}
#' 
#' @keywords internal
#' @noRd
rclust.specc <- function(input, k=2, kernel="rbfdot", kpar="automatic",
                         nystrom.red = FALSE, nystrom.sample = length(input$data)/6,
                         maxiter = 496, mod.sample = 0.75, na.action=na.omit,
                         type=c("extrinsic","intrinsic")){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if (is.matrix(input)){
    input = RiemBase::riemfactory(t(input), name="euclidean")
  }
  if ((class(input))!="riemdata"){
    stop("* rclust.specc : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  # acquire manifold name
  mfdname = tolower(input$name)
  # type of distance
  distype = match.arg(type)
  
  my.kernel = kernel
  my.kpar   = kpar
  my.ny.red = nystrom.red
  my.ny.sam = nystrom.sample
  my.iterations = maxiter
  my.mod.sample = mod.sample
  my.na.action  = na.omit
  
  #-------------------------------------------------------
  # equivariant embedding and apply "kernlab::specc"
  dataequiv  = rclust_equivariant(input) # now must be aligned in STATISTICAL convention
  kernoutput = kernlab::specc(dataequiv, k, kernel=my.kernel, kpar=my.kpar, 
                              nystrom.red=my.ny.red, nystrom.sample=my.ny.sam,
                              iterations=my.iterations, mod.sample=my.mod.sample,
                              na.action=my.na.action)
  
  # extract valid information to be reported later
  label      = as.integer(kernoutput)
  size       = kernoutput@size
  
  # borrow the internal function from "rstat.frechet.cube"
  frechet.internal = utils::getFromNamespace("rstat.frechet.cube","RiemStats")
  
  #-------------------------------------------------------
  # compute class-wise mean and BCSS 
  # (modify to correspond with WCSS: \sum |S_i| Var(S_i) <- |S_i| V_i)
  ulabel   = unique(label)
  datacube = rclust_stack3d(input)
  centers  = list()
  withinss = rep(0,length(ulabel))
  for (i in 1:length(ulabel)){
    idxi = which(label==ulabel[i]) # corresponding group elements
    tmpi = frechet.internal(rclust_cubesubset(datacube,idxi),mfdname,type=distype)
    centers[[i]] = tmpi$x
    withinss[i] = (length(idxi)*tmpi$variation)
  }
  wcss = sum(withinss)
  
  #-------------------------------------------------------
  # return output
  output = list()
  output$centers  = centers  # list of center matrices
  output$wcss     = wcss     # wcss as a value
  output$withinss = withinss # for each class
  output$label    = label    # class labels
  output$size     = size     # size of each cluster
  return(output)
}