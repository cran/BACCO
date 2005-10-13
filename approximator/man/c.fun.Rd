\name{c.fun}
\alias{c.fun}
\alias{cdash.fun}
\title{Correlations between points in parameter space}
\description{
Correlation matrices between (sets of) points in parameter space, both
prior (\code{c.fun()}) and posterior (\code{cdash.fun()}).
}
\usage{
c.fun(x, xdash=x, subsets, hpa)
cdash.fun(x, xdash=x, V=NULL, Vinv=NULL, D1, subsets, basis, hpa, method=2)
}
\arguments{
  \item{x,xdash}{Points in parameter space; or, if a matrix,
    interpret the rows as points in parameter space.  Note that the
    default value of \code{xdash} (viz \code{x}) will return the
    variance-covariance matrix of a set of points.}
  \item{D1}{Design matrix}
  \item{subsets}{Subset object}
  \item{hpa}{hyperparameter object}
  \item{basis}{Basis function}
  \item{V,Vinv}{In function \code{cdash.fun()}, the data covariance
    matrix and its inverse.  If
    \code{NULL}, the matrix will be calculated from scratch.  Supplying
    a precalculated value for \code{V}, and especially \code{Vinv},
    makes for very much faster execution (edepending on \code{method})}
  \item{method}{Integer specifying which of several algebraically
    identical methods to use.  See the source code for details, but
    default option 2 seems to be the best.  Bear in mind that option 3
    does not require inversion of a matrix, but is not faster in
    practice.}
}
\value{
Returns a matrix of covariances
}
\references{KOH2000}
\author{Robin K. S. Hankin}
\note{

  Do not confuse function \code{c.fun()}, which computes \eqn{c(x,x')}
  defined just below equation 7 on page 4 with \eqn{c_t(x,x')} defined
  in equation 3 on page 3.

  Consider the example given for two levels on page 4 just after
  equation 7:
  \eqn{c(x,x')=c_2(x,x')+\rho_1^2c_1(x,x')}{c(x,x')=c_2(x,x')+rho_1^2c_1(x,x')}
  is a kind of prior covariance matrix.  Matrix \eqn{c'(x,x')} is a
  posterior covariance matrix, conditional on the code observations.

  Function \code{A()} evaluates \eqn{c_t(x,x')} in a nice vectorized way.

  Equation 7 of KOH2000 contains a typo.
}
\seealso{\code{\link{A}}}
\examples{
data(toyapps)

x <- latin.hypercube(4,3)
rownames(x) <- c("ash" , "elm" , "oak", "pine")
xdash <- latin.hypercube(7,3)
rownames(xdash) <- c("cod","bream","skate","sole","eel","crab","squid")

cdash.fun(x=x,xdash=xdash, D1=D1.toy, basis=basis.toy,subsets=subsets.toy, hpa=hpa.toy)

# Now add a point whose top-level value is known:
x <- rbind(x,D1.toy[subsets.toy[[4]][1],])

cdash.fun(x=x,xdash=xdash, D1=D1.toy, basis=basis.toy,subsets=subsets.toy, hpa=hpa.toy)
# Observe how the bottom row is zero (up to rounding error)


}
\keyword{math} 