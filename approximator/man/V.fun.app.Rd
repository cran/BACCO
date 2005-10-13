\name{V.fun.app}
\alias{V.fun.app}
\title{Variance matrix}
\description{
Given a design matrix, a subsets object and a hyperparameter object,
return the variance matrix.  The \dQuote{app} of the function name means
\dQuote{approximator}, to distinguish it from function \code{V.fun()}
of the calibrator package.
}
\usage{
V.fun.app(D1, subsets, hpa)
}
\arguments{
  \item{D1}{Design matrix for level 1 code}
  \item{subsets}{Subsets object}
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
V.fun.app(D1.toy,subsets.toy,hpa.toy)
}
\keyword{array}
