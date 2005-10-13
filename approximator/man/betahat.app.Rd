\name{betahat.app}
\alias{betahat.app}
\alias{betahat.app.H}
\title{Estimate for beta}
\description{
Returns the estimator for beta; equation 5.  Function
\code{betahat.app()} returns the estimate in terms of fundamental
variables; \code{betahat.app.H()} requires the \code{H} matrix.
}
\usage{
betahat.app.H(H, V = NULL, Vinv = NULL, z)
betahat.app(D1, subsets, basis, hpa, z, use.Vinv=TRUE)
}
\arguments{
  \item{H}{In \code{betahat.app.H()}, the H matrix, eg that
    returned by \code{H.fun()}}
  \item{V}{Variance matrix}
  \item{Vinv}{Inverse of variance matrix.  If not supplied, it is
    calculated}
  \item{use.Vinv}{In function \code{betahat.app()}, a Boolean argument
    with default \code{TRUE} meaning to calculate the inverse of the
    \code{V} matrix; and \code{FALSE} meaning to use a method which does
    not involve calculating the inverse of \code{V}.  The default method
    seems to be faster; YMMV}
  \item{z}{vector of observations}
  \item{D1}{Design matrix for level 1 code}
  \item{subsets}{Subsets object}
  \item{basis}{Basis function}
  \item{hpa}{Hyperparameter object}
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)

 betahat.app(D1=D1.toy, subsets=subsets.toy, basis=basis.toy, hpa=hpa.toy, z=z.toy, use.Vinv=TRUE)

 H <- H.fun.app(D1=D1.toy, subsets=subsets.toy, basis=basis.toy,hpa=hpa.toy)
 V <- V.fun.app(D1=D1.toy, subsets=subsets.toy, hpa=hpa.toy)
 betahat.app.H(H=H,V=V,z=z.toy)

}
\keyword{array}
