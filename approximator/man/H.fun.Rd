\name{H.fun}
\alias{H.fun.app}
\title{The H matrix}
\description{
Returns the matrix of bases H.    The \dQuote{app} of the function name means
\dQuote{approximator}, to distinguish it from function \code{H.fun()}
of the calibrator package.
}
\usage{
H.fun.app(D1, subsets, basis, hpa)
}
\arguments{
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
H.fun.app(D1.toy , subsets=subsets.toy , basis=basis.toy , hpa=hpa.toy)
}
\keyword{array}
