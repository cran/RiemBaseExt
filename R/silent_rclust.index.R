#' Clustering Validation Index with Manifold-Valued Data
#' 
#' CVI
#' 
#' 
#' @keywords internal
#' @noRd
rclust.index <- function(input, label, 
                         index=c("CCC","CH","Cindex","DB","Dunn","Gamma","Gplus","McClain","Ptbiserial","Ratkowsky","SDbw","SDindex","Silhouette","Tau")){
  allindices = c("ch","silhouette","cindex","ptbiserial","db","mcclain","dunn","gamma","gplus","tau","ratkowsky","ccc","sdbw","sdindex")
  #-------------------------------------------------------
  if (is.matrix(input)){
    input = RiemBase::riemfactory(t(input), name="euclidean")
  }
  # must be of 'riemdata' class
  if ((class(input))!="riemdata"){
    stop("* rclust.index : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  if ((length(input$data)!=length(label))||(!is.vector(label))){
    stop("* rclust.index : length of a class label vector must correspond to that of provided data.")
  }
  # take care of label
  label = as.integer(label)
  # index type matching
  index = match.arg(tolower(index), allindices)
  
  
  #-------------------------------------------------------
  # computation
  # convert into data cube and extract mfdname
  dcube   = rclust_stack3d(input)
  mfdname = input$name
  memberlist = rclust_membership(label)
  if (index=="silhouette"){
    score = index_silhouette(dcube, mfdname, label, memberlist)
  } else if (index=="ch"){
    score = index_ch(dcube, mfdname, label, memberlist)
  } else if (index=="cindex"){
    score = index_cindex(dcube, mfdname, label, memberlist)
  } else if (index=="ptbiserial"){
    score = index_ptbiserial(dcube, mfdname, label, memberlist)
  } else if (index=="db"){
    score = index_db(dcube, mfdname, label, memberlist)
  } else if (index=="mcclain"){
    score = index_mcclain(dcube, mfdname, label, memberlist)
  } else if (index=="dunn"){
    score = index_dunn(dcube, mfdname, label, memberlist)
  } else if (index=="gamma"){
    score = index_gamma(dcube, mfdname, label, memberlist)
  } else if (index=="gplus"){
    score = index_gplus(dcube, mfdname, label, memberlist)
  } else if (index=="tau"){
    score = index_tau(dcube, mfdname, label, memberlist)
  } else if (index=="ratkowsky"){
    score = index_ratkowsky(dcube, mfdname, label, memberlist)
  } else if (index=="ccc"){
    score = index_ccc(dcube, mfdname, label, memberlist)
  } else if (index=="sdbw"){
    score = index_sdbw(dcube, mfdname, label, memberlist)
  } else if (index=="sdindex"){
    score = index_sdindex(dcube, mfdname, label, memberlist)
  }
  
  
  #-------------------------------------------------------
  # return score
  return(score)
}




############################################################################
# (01) Silhouette Index / Max / Rousseeux (1987)
#' @keywords internal
#' @noRd
index_silhouette <- function(dcube, mfdname, label, memberlist){
  # compute 'pdist'
  distmat = rclust_pdist_cube(dcube, mfdname, type="intrinsic")
  # labeling care
  ulabel = unique(label)
  k      = length(ulabel)
  if (k==1){
    stop("* rclust.index : for a single-label clustering, Silhouette index is not defined.")
  }
  # let's iterate, take a lot of time !
  n = length(label)
  vec.a = rep(0,n)
  vec.b = rep(0,n)
  for (i in 1:n){ # for each data object
    # label of i-th element
    idlabel = which(ulabel==label[i])
    # compute a(i)
    idx.samelabel = setdiff(memberlist[[idlabel]],i)
    vec.a[i] = base::mean(as.vector(distmat[i,idx.samelabel]))
    # compute b(j)
    otherlabels = setdiff(ulabel, label[i])
    dics = rep(0,k-1)
    for (j in 1:(k-1)){
      jlabel = which(ulabel==otherlabels[j])
      idx.jlabel = memberlist[[jlabel]]
      dics[j] = base::mean(as.vector(distmat[i,idx.jlabel]))
    }
    vec.b[i] = min(dics)
  }
  # now compute using vectorized operations
  vec.s = ((vec.b-vec.a)/(base::pmax(vec.a, vec.b)))
  # return
  return(mean(vec.s))
}

############################################################################
# (02) CH / Max / Calinski and Harabasz (1974)
#' @keywords internal
#' @noRd
index_ch <- function(dcube, mfdname, label, memberlist){
  # basic parameters
  n = length(label)
  q = length(unique(label))
  if (q==1){
    stop("* rclust.index : for a single-label clustering, CH index is not defined.")
  }
  
  ## Previous Computation : Locally at the Tangent Space
  # Bq = rclust_index_Bq(dcube, label, mfdname, memberlist)
  # Wq = rclust_index_Wq(dcube, label, mfdname, memberlist)
  # term.num = sum(diag(Bq))/(q-1)
  # term.den = sum(diag(Wq))/(n-1)
  
  # # compute preliminary ones
  BGSS = rclust_index_BGSS(dcube, label, mfdname, memberlist)
  WGSS = rclust_index_WGSS(dcube, label, mfdname, memberlist)
  # compute the score
  term.num = BGSS/(q-1)
  term.den = WGSS/(n-q)
  score    = term.num/term.den
  # return
  return(score)
}

############################################################################
# (03) Cindex / Min / Hubert and Levin (1976)
#' @keywords internal
#' @noRd
index_cindex <- function(dcube, mfdname, label, memberlist){
  # basic parameters
  n  = length(label)
  q  = length(unique(label))
  
  Ns = rclust_index_Ns(label, memberlist) # Nt, Ns, Nw at once
  Nw = Ns$Nw
  
  # compute pairwise distance matrix
  dmat = rclust_pdist_cube(dcube, mfdname, type="intrinsic")
  # compute Sw
  Sw = rclust_index_Sw(dmat, label, mfdname, memberlist)
  # Smin and Smax
  vecd = sort(as.vector(dmat[lower.tri(dmat)])) # sort the distances in an increasing order
  Smin = sum(head(vecd, Nw))
  Smax = sum(tail(vecd, Nw))
  # compute score and return output
  Cindex = ((Sw-Smin)/(Smax-Smin))
  return(Cindex)
}


############################################################################
# (04) Ptbiserial / Max / Milligan 1980
#' @keywords internal
#' @noRd
index_ptbiserial <- function(dcube, mfdname, label, memberlist){
  # basic parameters
  n = length(label)
  # info : counts
  Ns = rclust_index_Ns(label, memberlist)
  Nb = as.double(Ns$Nb)
  Nw = as.double(Ns$Nw)
  Nt = as.double(Ns$Nt)
  # info : Sw- and Sb-corresponding bar matrices
  dmat   = rclust_pdist_cube(dcube, mfdname, type="intrinsic")
  Sw.bar = rclust_index_Sw(dmat, label, mfdname, memberlist)/Nw
  Sb.bar = rclust_index_Sb(dmat, label, mfdname, memberlist)/Nb
  # compute
  Sd    = stats::sd(as.vector(dmat[lower.tri(dmat)])) # standard deviation
  score = (((Sb.bar-Sw.bar)*sqrt(Nw*Nb/(Nt^2)))/Sd)
  return(score)
}

############################################################################
# (05) DB / Min / Daviews and Bouldin 1979
#' @keywords internal
#' @noRd
index_db <- function(dcube, mfdname, label, memberlist){
  # label 
  ulabel = unique(label)
  q      = length(ulabel)
  dnrow  = nrow(dcube)
  dncol  = ncol(dcube)
  # compute class-wise mean (cube)
  mean.class = rclust_classmean(dcube, label, mfdname, memberlist)
  # compute distance between centroids
  dmat.centroids = rclust_pdist_cube(mean.class, mfdname, type="intrinsic")
  # compute in-class distance measures
  dvec.inclass = rep(0,q)
  for (i in 1:q){
    idlabel = memberlist[[which(ulabel==ulabel[i])]]
    dcube1  = array(0,c(dnrow,dncol,1)); dcube1[,,1] = mean.class[,,i]
    dcube2  = rclust_cubesubset(dcube, idlabel)
    dvec.inclass[i] = stats::sd(as.vector(rclust_pdist2_cube(dcube1, dcube2, mfdname, type="intrinsic")))
  }
  # now, compute the score
  vec.db = rep(0,q)
  for (k in 1:q){
    idrest = setdiff(1:q,k)
    del.k = dvec.inclass[k]
    del.l = dvec.inclass[-k] # vector
    dvecs = as.vector(dmat.centroids[k,idrest])
    vec.db[k] = max((rep(del.k,(q-1))+del.l)/dvecs)
  }
  return(mean(vec.db))
}

############################################################################
# (07) McClain / Min / McClain and Rao 1975
#' @keywords internal
#' @noRd
index_mcclain <- function(dcube, mfdname, label, memberlist){
  # info : counts
  Ns = rclust_index_Ns(label, memberlist)
  Nb = as.double(Ns$Nb)
  Nw = as.double(Ns$Nw)
  Nt = as.double(Ns$Nt)
  # info : Sw- and Sb-corresponding bar matrices
  dmat   = rclust_pdist_cube(dcube, mfdname, type="intrinsic")
  Sw.bar = rclust_index_Sw(dmat, label, mfdname, memberlist)/Nw
  Sb.bar = rclust_index_Sb(dmat, label, mfdname, memberlist)/Nb
  # output
  output = Sw.bar/Sb.bar
  return(output)
}

############################################################################
# (08) Dunn / Max / Dunn 1974
#' @keywords internal
#' @noRd
index_dunn <- function(dcube, mfdname, label, memberlist){
  # basic information
  q      = length(memberlist) # number of unique clusters
  dmat   = rclust_pdist_cube(dcube, mfdname, type="intrinsic") # distance matrix
  ulabel = unique(label)
  # compute 1 : diameter
  vec.diam = rep(0,q)
  for (i in 1:q){
    tgt         = memberlist[[which(ulabel==ulabel[i])]]
    vec.diam[i] = max(as.vector(dmat[tgt,tgt]))
  }
  # compute 2 : inter-cluster distance
  mat.dist = array(0,c(q,q))
  for (i in 1:(q-1)){
    id1 = memberlist[[which(ulabel==ulabel[i])]]
    for (j in (i+1):q){
      id2   = memberlist[[which(ulabel==ulabel[j])]]
      mat.dist[i,j] = min(as.vector(dmat[id1,id2]))
    }
  }
  
  # compute : final step
  term.den = max(vec.diam)
  term.num = min(setdiff(unique(as.vector(mat.dist)),0))
  output   = (term.num/term.den)
  return(output)
  
}

############################################################################
# (09) Gamma / Max / Baker and Hubert 1975
#' @keywords internal
#' @noRd
index_gamma <- function(dcube, mfdname, label, memberlist){
  # compute distance matrix
  dmat   = rclust_pdist_cube(dcube, mfdname, type="intrinsic") # distance matrix
  # compute concordance information
  concordance = rclust_concordant(dmat, label, memberlist)
  sp = as.double(concordance$con)
  sm = as.double(concordance$dis)
  # score
  output = (sp-sm)/(sp+sm)
  return(output)
}

############################################################################
# (10) Gplus / Min / Rohlf 1974
#' @keywords internal
#' @noRd
index_gplus <- function(dcube, mfdname, label, memberlist){
  # compute distance matrix
  dmat   = rclust_pdist_cube(dcube, mfdname, type="intrinsic") # distance matrix
  # compute concordance information
  concordance = rclust_concordant(dmat, label, memberlist)
  sm = as.double(concordance$dis)
  # compute Nt
  n  = length(label)
  Nt = as.double(n*(n-1)/2)
  # compute score
  output = (2*sm/(Nt*(Nt-1)))
  return(output)
}

############################################################################
# (11) Tau / Max / Rohlf 1974
#' @keywords internal
#' @noRd
index_tau <- function(dcube, mfdname, label, memberlist){
  # compute distance matrix
  dmat   = rclust_pdist_cube(dcube, mfdname, type="intrinsic") # distance matrix
  # compute concordance information
  concordance = rclust_concordant(dmat, label, memberlist)
  sp = as.double(concordance$con)
  sm = as.double(concordance$dis)
  # compute Nt, Nb, Nw
  info.n = rclust_index_Ns(label, memberlist)
  Nt = as.double(info.n$Nt)
  Nb = as.double(info.n$Nb)
  Nw = as.double(info.n$Nw)
  
  # compute score
  term.num = (sp-sm)
  term.den = sqrt(Nb*Nw*Nt*(Nt-1)/2)
  output   = term.num/term.den
  return(output)
}

############################################################################
# (12) Ratkowsky / Max / Ratkowsky and Lance 1978
#' @keywords internal
#' @noRd
index_ratkowsky <- function(dcube, mfdname, label, memberlist){
  # ulabel
  ulabel = unique(label)
  q      = length(memberlist)
  # equivariant embedding first.
  X = rclust_equivariant_cube(dcube, mfdname)
  p = ncol(X)
  # compute n_k vector
  vec.Nk = rep(0,q)
  for (i in 1:q){
    vec.Nk[i] = length(memberlist[[which(ulabel==ulabel[i])]])
  }
  # compute class-wise and global mean
  mean.class = array(0,c(q,p))
  for (i in 1:q){
    mean.class[i,] = as.vector(colMeans(X[memberlist[[which(ulabel==ulabel[i])]],]))
  }
  mean.global = as.vector(colMeans(X))
  # compute 1 : BGSS
  BGSS = rep(0,p)
  for (j in 1:p){
    BGSS[j] = sum(((as.vector(mean.class[,j])-mean.global[j])^2)*vec.Nk)
  }
  # compute 2 : TSS
  TSS = rep(0,p)
  for (i in 1:p){
    TSS[i] = sum((as.vector(X[,i])-mean.global[i])^2)
  }
  # compute the scoer
  return((sqrt((1/p)*sum(BGSS/TSS)))/sqrt(q))
}

############################################################################
# (13) CCC / Max / Sarle 1983
#' @keywords internal
#' @noRd
index_ccc <- function(dcube, mfdname, label, memberlist){
  # parameters and equivariant embedding
  ulabel = unique(label)
  X      = rclust_equivariant_cube(dcube, mfdname) # extrinsic !
  n      = nrow(X)
  p      = ncol(X)
  q      = length(memberlist)
  
  # preliminary computation
  XtX    = t(X)%*%X
  vec.s  = eigen(XtX/(n-1))$values      # squared-root eigenvalues from empirical gram matrix
  print(vec.s)
  vec.s[(vec.s<=0)]=sqrt(.Machine$double.eps); vec.s = sqrt(vec.s);
  val.v = prod(vec.s)
  val.c = (val.v/q)^(1/p)
  vec.u = (vec.s/val.c)
  val.p = max(pmin(which(vec.u>=1),(q-1)))
  Z     = rclust_Z(label, memberlist)
  X.bar = base::solve((t(Z)%*%Z),(t(Z)%*%X))
  
  # main part 1 : R2 and ER2
  R2  = 1-(sum(diag(XtX - (t(X.bar)%*%t(Z)%*%Z%*%X.bar)))/sum(diag(XtX)))
  ER2 = 1-((((sum(1/(n+vec.u[1:val.p]))) + (sum((as.vector(vec.u[(val.p+1):p])^2)/(n+vec.u[(val.p+1):p]))))/(sum(vec.u^2)))*(((n-q)^2)/n)*(1+(4/n)))
  
  # main part 2 : truly main
  output = log((1-ER2)/(1-R2))*sqrt(n*val.p/2)/((0.001+ER2)^(1.2))
  return(output)
}

############################################################################
# (14) SDindex / Min / Halkidi et al 2000
#' @keywords internal
#' @noRd
index_sdindex <- function(dcube, mfdname, label, memberlist){
  # parameters and equivariant embedding
  ulabel = unique(label)
  X      = rclust_equivariant_cube(dcube, mfdname) # extrinsic !
  q      = length(memberlist)
  
  # compute 1 : Scat(q)
  sig.local = list()
  for (i in 1:q){
    sig.local[[i]] = apply(X[memberlist[[which(ulabel==ulabel[i])]],], 2, var)
  }
  sig.global = apply(X, 2, var)
  Scat = 0
  for (i in 1:q){
    tgt  = sig.local[[i]]
    Scat = Scat + sqrt(sum(tgt^2))
  }
  Scat = Scat/(sqrt(sum(as.vector(sig.global)^2))*q)
  
  # compute 2 : Dis
  classmean = rclust_classmean(dcube, label, mfdname, memberlist)
  distmat   = rclust_pdist_cube(classmean, mfdname, type="intrinsic")
  halfD     = as.vector(distmat[upper.tri(distmat)])
  Dmax      = max(halfD)
  Dmin      = min(halfD)
  
  Dis = 0.0
  for (k in 1:q){
    Dis = Dis + (1/sum(as.vector(distmat[k,])))
  }
  Dis = Dis*Dmax/Dmin
  
  # output
  output = Dis*(1+Scat)
  return(output)
}

############################################################################
# (15) SDbw / Min / Halkidi and Vazirgiannis 2001
#' @keywords internal
#' @noRd
index_sdbw <- function(dcube, mfdname, label, memberlist){
  # parameters and equivariant embedding
  ulabel = unique(label)
  X      = rclust_equivariant_cube(dcube, mfdname) # extrinsic !
  q      = length(memberlist)
  
  # compute 1 : Scat(q)
  sig.local = list()
  for (i in 1:q){
    sig.local[[i]] = apply(X[memberlist[[which(ulabel==ulabel[i])]],], 2, var)
  }
  sig.global = apply(X, 2, var)
  Scat = 0
  for (i in 1:q){
    tgt  = sig.local[[i]]
    Scat = Scat + sqrt(sum(tgt^2))
  }
  Scat = Scat/(sqrt(sum(as.vector(sig.global)^2))*q)
  
  # compute 2 : Density.bw with Stdev : refer to clusterCrit
  # 2-1. threshold
  Stdev = 0
  for (k in 1:q){
    Stdev = Stdev + sqrt(sum(sig.local[[k]]^2))
  }
  Stdev = sqrt(Stdev)/q
  # 2-2. classwise mean in extrinsic manner
  classmean = rclust_equivariant_cube(rclust_classmean(dcube, label, mfdname, memberlist), mfdname)
  # 2-3. iterate
  density.bw = 0
  for (i in 1:(q-1)){
    id1 = memberlist[[which(ulabel==ulabel[i])]]
    for (j in (i+1):q){
      id2 = memberlist[[which(ulabel==ulabel[j])]]
      
      idunion = c(id1,id2) # join two groups
      cluster1 = as.vector(classmean[i,])
      cluster2 = as.vector(classmean[j,])
      centroid = as.vector((cluster1+cluster2)/2)
      
      val.top  = rclust_index_density(centroid, X[idunion,], Stdev)
      val.bot1 = rclust_index_density(cluster1, X[idunion,], Stdev)
      val.bot2 = rclust_index_density(cluster2, X[idunion,], Stdev)
      if (max(c(val.bot1,val.bot2))>0){
        density.bw = density.bw + (val.top/max(c(val.bot1,val.bot2)))
      }
    }
  }
  density.bw = density.bw*2/(q*(q-1))
  
  # return output
  return(Scat+density.bw)
}

############################################################################
############################################################################
## Et Cetera
# * Scott / Maximal Incremental / Scott and Symons 1971
#' @keywords internal
#' @noRd
index_scott <- function(dcube, mfdname, label, memberlist){
  # compute Wq
  Wq = rclust_index_Wq(dcube, label, mfdname)
  # comptue T (as Tmat)
  meanfunc = utils::getFromNamespace("rbase.mean.cube","RiemBase")
  mean.global = meanfunc(dcube, mfdname)$x
  log.pulled = cpp_dispersion(dcube, mean.global, mfdname)
  Tmat       = (log.pulled%*%t(log.pulled))
  # compute the score
  n = length(label)
  score = log(det(Tmat)/det(Wq))*n
  return(score)
}

# for IRIS dataset [,2:4], McClain and CCC are reporting weird results.