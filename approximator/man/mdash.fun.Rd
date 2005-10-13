\name{mdash.fun}
\alias{mdash.fun}
\title{Mean of Gaussian process}
\description{
Returns the mean of the Gaussian process conditional on the observations
and the hyperparameters
}
\usage{
mdash.fun(x, D1, subsets, hpa, Vinv = NULL, use.Vinv = TRUE, z, basis)
}

\arguments{
  \item{x}{Point at which mean is desired}
  \item{D1}{Code design matrix for level 1 code}
  \item{subsets}{subsets object}
  \item{hpa}{Hyperparameter object}
  \item{Vinv}{Inverse of the variance matrix; if \code{NULL}, the
  function will calculate it}
  \item{use.Vinv}{Boolean, with default \code{TRUE} meaning to use the
  inverse of V and \code{FALSE} meaning to use a method that does not
  involve inverting V}
  \item{z}{observations}
  \item{basis}{Basis functions}
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)
mdash.fun(x=1:3,D1=D1.toy,subsets=subsets.toy,hpa=hpa.toy,z=z.toy,basis=basis.toy)

uu <- rbind(1:3,1,3:1,1:3)
rownames(uu) <- c("first","second","third","fourth")

mdash.fun(x=uu,D1=D1.toy,subsets=subsets.toy,hpa=hpa.toy,z=z.toy,basis=basis.toy)

}
\keyword{array}
