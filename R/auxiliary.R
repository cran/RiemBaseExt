## auxiliary functions
#   (1) aux_stack3d : given 'object.riem' class, stack the data as 3d
#   (2) aux_eqembed : equivariant embedding (stacked as rows for future use)


# (1) aux_stack3d ---------------------------------------------------------
#' @keywords internal
#' @noRd
aux_stack3d <- function(object.riem){
  msize = object.riem$size
  ndata = length(object.riem$data)
  
  matdata = array(0,c(msize[1], msize[2], ndata))
  for (i in 1:ndata){
    tgt = object.riem$data[[i]]
    if (is.vector(tgt)){
      matdata[,,i] = as.matrix(tgt)
    } else {
      matdata[,,i] = tgt
    }
  }
  return(matdata)
}


# (2) aux_eqembed ---------------------------------------------------------
#' @keywords internal
#' @noRd
aux_eqembed <- function(object.riem){
  datacube = aux_stack3d(object.riem)
  datamat  = t(convert_3toM(datacube, object.riem$name))
  return(datamat)
}