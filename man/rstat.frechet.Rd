% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstat.frechet.R
\name{rstat.frechet}
\alias{rstat.frechet}
\title{Fréchet Mean and Variation of Manifold-valued Data}
\usage{
rstat.frechet(
  input,
  type = c("intrinsic", "extrinsic"),
  int.eps = 1e-06,
  parallel = FALSE
)
}
\arguments{
\item{input}{a S3 object of \code{riemdata} class. See \code{\link[RiemBase]{riemfactory}} for more details.}

\item{type}{type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.}

\item{int.eps}{stopping criterion for the norm of gradient.}

\item{parallel}{a flag for enabling parallel computation.}
}
\value{
a named list containing
\describe{
\item{mu}{an estimated Fréchet mean matrix.}
\item{variation}{Fréchet variation with the estimated mean.}
}
}
\description{
For manifold-valued data, Fréchet mean is the solution of following cost function,
\deqn{\textrm{min}_x \sum_{i=1}^n \rho^2 (x, x_i),\quad x\in\mathcal{M}}
for a given data \eqn{\{x_i\}_{i=1}^n} and \eqn{\rho(x,y)} is the geodesic distance 
between two points on manifold \eqn{\mathcal{M}}. It uses a gradient descent method 
with a backtracking search rule for updating.
}
\examples{
### Generate 50 data points on Sphere S^2 near (0,0,1).
ndata = 50
theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
tmpy  = sin(theta) + rnorm(ndata,sd=0.1)

### Wrap it as 'riemdata' class
data  = list()
for (i in 1:ndata){
  tgt = c(tmpx[i],tmpy[i],1)
  data[[i]] = tgt/sqrt(sum(tgt^2)) # project onto Sphere
}
data = RiemBase::riemfactory(data, name="sphere")

### Compute Fréchet Mean and Variation
out1 = rstat.frechet(data)                   # intrinsic
out2 = rstat.frechet(data,parallel=TRUE)     # parallel implementation
out3 = rstat.frechet(data, type="extrinsic") # extrinsic

}
\references{
\insertRef{karcher_riemannian_1977}{RiemBaseExt}

\insertRef{kendall_probability_1990}{RiemBaseExt}

\insertRef{afsari_convergence_2013}{RiemBaseExt}
}
