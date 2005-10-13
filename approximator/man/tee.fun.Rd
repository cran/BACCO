\name{tee.fun}
\alias{tee.fun}
\title{Returns generalized distances}
\description{
Returns generalized distances from a point to the design matrix as per
equation 10
}
\usage{
tee.fun(x, D1, subsets, hpa)
}
\arguments{
  \item{x}{Point in parameter space}
  \item{D1}{Design matrix for level 1 code}
  \item{subsets}{subsets object}
  \item{hpa}{Hyperparameter object}
}
\details{
See equation 10
}

\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{

data(toyapps)
tee.fun(x=1:3, D1=D1.toy, subsets=subsets.toy, hpa=hpa.toy)
}
\keyword{array}
