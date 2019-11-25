############# Auxiliary R Functions #########################
# (00) rclust_mat2cube      : transform a matrix into 3d cube
# (01) rclust_stack3d       : stack data into 3d cube
#      rclust_stack3d2      : stack two data into one 3d cube
# (02) rclust_equivariant   : make an equivariant embedding stacked as rows
#      rclust_equivariant_cube
# (03) rclust_pdist         : pairwise intrinsic/extrinsic distances
#      rclust_pdist_cube  
# (04) rclust_pdist2        : pairwise intrinsic/extrinsic distances for two-set data
#      rclust_pdist2_cube
# (05) rclust_cubesubset    : remove anomaly in subsetting cube by third dimension
# (06) rclust_classmean     : compute class-wise mean cube
# (07) rclust_index_Bq      : between-group                    : CH is not okay. Change to Distance Measures
#      rclust_index_Wq      : within -group dispersion matrix
#      rclust_index_WGSS    : within -group sum of squared distances
#      rclust_index_BGSS    : between-group 
# (08) rclust_index_Sw      : within-cluster
#      rculst_index_Sb      : between-cluster sum of distances given a distance matrix
# (09) rclust_index_Ns      : number of counts information
# (10) rclust_membership    : membership of clusterings
# (11) rclust_concordant    : s(+) and s(-)
# (12) rclust_Z             : membership matrix of size (n x q)
# (13) rclust_index_density : use for SDbw method

# (00) rclust_mat2cube    : transform a matrix into 3d cube
#' @keywords internal
#' @noRd
rclust_mat2cube <- function(inputmat){
  output = array(0,c(nrow(inputmat),ncol(inputmat),1))
  output[,,1] = inputmat
  return(output)
}

# (01) rclust_stack3d     : stack data into 3d cube
#' @keywords internal
#' @noRd
rclust_stack3d <- function(riemdata){
  msize = riemdata$size
  ndata = length(riemdata$data)
  
  matdata = array(0,c(msize[1], msize[2], ndata))
  for (i in 1:ndata){
    matdata[,,i] = (riemdata$data[[i]])
  }
  return(matdata)
}
#' @keywords internal
#' @noRd
rclust_stack3d2 <- function(riem1, riem2){
  msize = riem1$size
  ndat1 = length(riem1$data)
  ndat2 = length(riem2$data)
  
  matdata = array(0,c(msize[1], msize[2], (ndat1+ndat2)))
  for (i in 1:ndat1){
    matdata[,,i] = (riem1$data[[i]])
  }
  for (i in (ndat1+1):(ndat1+ndat2)){
    matdata[,,i] = (riem2$data[[i-ndat1]])
  }
  return(matdata)
}

# (02) rclust_equivariant : make an equivariant embedding stacked as rows
#' @keywords internal
#' @noRd
rclust_equivariant <- function(input){
  mfddata = rclust_stack3d(input)
  mfdname = input$name
  output  = cpp_equivariant(mfddata, mfdname)
  return(output)
}
#' @keywords internal
#' @noRd
rclust_equivariant_cube <- function(dcube, mfdname){
  return(cpp_equivariant(dcube, mfdname))
}

# (03) rclust_pdist       : pairwise intrinsic/extrinsic distances
#' @keywords internal
#' @noRd
rclust_pdist <- function(input, type=c("extrinsic","intrinsic")){
  mfddata = rclust_stack3d(input)
  mfdname = input$name
  distype = match.arg(type)
  
  if (distype=="intrinsic"){
    output = cpp_distint(mfddata, mfdname)
  } else {
    output = cpp_distext(mfddata, mfdname)
  }
  return(output)
}
#' @keywords internal
#' @noRd
rclust_pdist_cube <- function(datacube, mfdname, type=c("extrinsic","intrinsic")){
  mfddata = datacube
  distype = match.arg(type)
  
  if (distype=="intrinsic"){
    output = cpp_distint(mfddata, mfdname)
  } else {
    output = cpp_distext(mfddata, mfdname)
  }
  return(output)
}
# (04) rclust_pdist2      : pairwise intrinsic/extrinsic distances for two-set data
#' @keywords internal
#' @noRd
rclust_pdist2 <- function(input, type=c("extrinsic","intrinsic")){
  mfddata = rclust_stack3d(input)
  mfdname = input$name
  distype = match.arg(type)
  
  if (distype=="intrinsic"){
    output = cpp_dist2int(mfddata, mfdname)
  } else {
    output = cpp_dist2ext(mfddata, mfdname)
  }
  return(output)
}
#' @keywords internal
#' @noRd
rclust_pdist2_cube <- function(cube1, cube2, mfdname, type=c("extrinsic","intrinsic")){
  distype = match.arg(type)
  if (distype=="intrinsic"){
    output = cpp_dist2int(cube1, cube2, mfdname)
  } else {
    output = cpp_dist2ext(cube1, cube2, mfdname)
  }
  return(output)
}

# (05) rclust_cubesubset  : remove anomaly in subsetting cube by third dimension
#' @keywords internal
#' @noRd
rclust_cubesubset <- function(datacube, subsetid){
  nrow = dim(datacube)[1]
  ncol = dim(datacube)[2]
  nslice = length(subsetid)
  
  output = array(0,c(nrow,ncol,nslice))
  for (i in 1:nslice){
    tmp = datacube[,,subsetid[i]]
    if (is.vector(tmp)){
      tmp = as.matrix(tmp)
    }
    output[,,i]=tmp
  }
  return(output)
}

# (06) rclust_classmean   : compute class-wise mean cube
#' @keywords internal
#' @noRd
rclust_classmean <- function(datacube, label, mfdname, memberlist){
  # get RiemBase::rbase.mean.cube
  meanfunc = utils::getFromNamespace("rbase.mean.cube","RiemBase")
  ulabel = unique(label)
  classmean = array(0,c(nrow(datacube),ncol(datacube),length(ulabel)))
  for (i in 1:length(ulabel)){
    idx.i = memberlist[[which(ulabel==ulabel[i])]]
    if (length(idx.i)==1){
      classmean[,,i] = rclust_cubesubset(datacube, idx.i)
    } else {
      classmean[,,i] = meanfunc(rclust_cubesubset(datacube, idx.i), mfdname)$x
    }
  }
  return(classmean)
}

# (07) rclust_index_Bq    : between-group
#      rclust_index_Wq    : within -group dispersion matrix
#' @keywords internal
#' @noRd
rclust_index_Bq <- function(datacube, label, mfdname, memberlist){
  # unique label
  ulabel = unique(label)
  # get RiemBase::rbase.mean.cube
  meanfunc = utils::getFromNamespace("rbase.mean.cube","RiemBase")
  # compute intrinsic mean for global and local/per-class ones
  mean.global = meanfunc(datacube, mfdname)$x
  mean.local  = rclust_classmean(datacube, label, mfdname, memberlist)
  # count the size of clusters
  vec.count   = rep(0,length(ulabel))
  for (i in 1:length(ulabel)){
    vec.count[i] = length(memberlist[[which(ulabel==ulabel[i])]])
  }
  # compute log-pulled vectors of class means : stacked as columns
  log.pulled = cpp_dispersion(mean.local, mean.global, mfdname)
  dimsize    = nrow(log.pulled)
  # compute Bq
  Bq = array(0,c(dimsize,dimsize))
  for (i in 1:ncol(log.pulled)){
    Bq = Bq+((base::outer(log.pulled[,i],log.pulled[,i]))*vec.count[i])
  }
  # return output
  return(Bq)
}
#' @keywords internal
#' @noRd
rclust_index_BGSS <- function(datacube, label, mfdname, memberlist){
  ulabel = unique(label)
  meanfunc = utils::getFromNamespace("rbase.mean.cube","RiemBase")
  # means of clusters and global one
  mean.global = rclust_mat2cube(meanfunc(datacube, mfdname)$x)
  mean.local  = rclust_classmean(datacube, label, mfdname, memberlist)
  # compute distances from cluster means to global mean
  vec.dist    = as.vector(cpp_dist2int(mean.global, mean.local, mfdname))
  vec.count   = rep(0,length(memberlist))
  for (i in 1:length(memberlist)){
    vec.count[i] = length(memberlist[[which(ulabel==ulabel[i])]])
  }
  # compute output
  output = sum((vec.dist^2)*vec.count)
  return(output)
}
#' @keywords internal
#' @noRd
rclust_index_Wq <- function(datacube, label, mfdname, memberlist){
  ulabel = unique(label)
  # get RiemBase::rbase.mean.cube
  meanfunc = utils::getFromNamespace("rbase.mean.cube","RiemBase")
  # STEP 1 : initial one should be done by myself
  #   label selection
  idlabel    = memberlist[[which(ulabel==ulabel[1])]]
  class.mean = meanfunc(rclust_cubesubset(datacube,idlabel), mfdname)$x
  class.mems = rclust_cubesubset(datacube,idlabel)
  log.pulled = cpp_dispersion(class.mems, class.mean, mfdname)
  Wq         = (log.pulled%*%t(log.pulled))
  
  # STEP 2 : rest classes  
  if (length(ulabel)!=1){
    for (i in 2:length(ulabel)){
      idlabel    = memberlist[[which(ulabel==ulabel[i])]]
      class.mean = meanfunc(rclust_cubesubset(datacube,idlabel), mfdname)$x
      class.mems = rclust_cubesubset(datacube,idlabel)
      log.pulled = cpp_dispersion(class.mems, class.mean, mfdname)
      Wq         = Wq + (log.pulled%*%t(log.pulled))
    }
  } 
  # return output
  return(Wq)
}
#' @keywords internal
#' @noRd
rclust_index_WGSS <- function(datacube, label, mfdname, memberlist){
  ulabel = unique(label)
  q      = length(memberlist)
  # get RiemBase::rbase.mean.cube
  meanfunc = utils::getFromNamespace("rbase.mean.cube","RiemBase")
  #   label selection
  distvec    = rep(0,q)
  for (i in 1:q){
    idlabel = memberlist[[which(ulabel==ulabel[i])]]
    class.mean = rclust_mat2cube(meanfunc(rclust_cubesubset(datacube,idlabel), mfdname)$x)
    class.mems = rclust_cubesubset(datacube,idlabel)
    distvec[i] = sum(as.vector(cpp_dist2int(class.mean, class.mems, mfdname))^2)
  }
  return(sum(distvec))
}


# (08) sum of distances ---------------------------------------------------
#   rclust_index_Sw : within-cluster
#   rclust_index_Sb : between-cluster sum of distances given a distance matrix
#' @keywords internal
#' @noRd
rclust_index_Sw <- function(distmat, label, mfdname, memberlist){
  ulabel = unique(label)
  q      = length(ulabel)
  output = 0
  for (i in 1:q){
    idlabel = memberlist[[which(ulabel==ulabel[i])]]
    if (length(idlabel)!=1){
      partmat = distmat[idlabel,idlabel]
      output  = output + sum(as.vector(partmat[lower.tri(partmat)]))
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
rclust_index_Sb <- function(distmat, label, mfdname, memberlist){
  ulabel = unique(label)
  q      = length(ulabel)
  score  = 0.0
  for (k in 1:(q-1)){
    label.k = memberlist[[which(ulabel==ulabel[k])]]
    for (j in (k+1):q){
      label.j = memberlist[[which(ulabel==ulabel[j])]]
      score   = score + sum(as.vector(distmat[label.k, label.j]))
    }
  }
  return(score)
}

# (09) rclust_index_Ns    : number of counts information
#' @keywords internal
#' @noRd
rclust_index_Ns <- function(label, memberlist){
  # basic parameters
  ulabel = unique(label)
  n = length(label)
  q = length(ulabel)
  # Nt
  Nt = as.integer(n*(n-1)/2)
  # Nw
  Nw = 0
  for (i in 1:q){
    Nk = length(memberlist[[which(ulabel==ulabel[i])]])
    Nw = as.integer(Nk*(Nk-1)/2) + Nw
  }
  # Nb
  Nb = Nt - Nw
  # return output
  output = list()
  output$Nt = Nt
  output$Nb = Nb
  output$Nw = Nw
  return(output)
}

# (10) rclust_membership  : membership of clusterings
#' @keywords internal
#' @noRd
rclust_membership <- function(label){
  ulabel = unique(label)
  memvec = list()
  for (i in 1:length(ulabel)){
    memvec[[i]] = which(label==ulabel[i])
  }
  return(memvec)
}

# (11) rclust_concordant  : s(+) and s(-)
#' @keywords internal
#' @noRd
rclust_concordant <- function(distmat, label, membership){
  # get some parameters
  q      = length(membership)
  ulabel = unique(label)
  # 1. compute all within-distance vector
  vec.within = c()
  for (i in 1:q){
    idx = membership[[which(ulabel==ulabel[i])]]
    partmat = distmat[idx,idx]
    vec.within = c(vec.within, as.vector(partmat[upper.tri(partmat)]))
  }
  # 2. compute all between-distance vector
  vec.between = c()
  for (i in 1:(q-1)){
    id1 = membership[[which(ulabel==ulabel[i])]]
    for (j in (i+1):q){
      id2 = membership[[which(ulabel==ulabel[j])]]
      partmat = distmat[id1,id2]
      vec.between = c(vec.between, as.vector(partmat))
    }
  }
  # 3. count all numbers
  niter = length(vec.within)
  count.con = 0
  count.dis = 0
  for (i in 1:niter){
    count.con = count.con + sum((vec.between>vec.within[i]))
    count.dis = count.dis + sum((vec.between<vec.within[i]))
  }
  # return
  output = list()
  output$con = count.con
  output$dis = count.dis
  return(output)
}

# (12) rclust_Z           : membership matrix of size (n x q)
#' @keywords internal
#' @noRd
rclust_Z <- function(label, membership){
  n = length(label)
  q = length(membership)
  ulabel = unique(label)
  
  output = array(0,c(n,q))
  for (i in 1:q){
    id = membership[[which(ulabel==ulabel[i])]]
    output[id,i] = 1
  }
  return(output)
}

# (13) rclust_index_density : use for SDbw method in extrinsic manner
#' @keywords internal
#' @noRd
rclust_index_density <- function(meanvec, datamat, thr){
  centered = datamat - matrix(rep(meanvec,nrow(datamat)),nrow=nrow(datamat),byrow=TRUE)
  distvec  = apply(centered, 1, function(x){sqrt(sum(x^2))})
  return(sum(distvec<thr))
}