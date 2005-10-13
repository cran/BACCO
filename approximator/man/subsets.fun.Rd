\name{subsets.fun}
\alias{is.nested}
\alias{is.strict}
\alias{subsets.fun}
\title{Generate and test subsets}
\description{
  Create a list of subsets (\code{subsets.fun()}); or, 
given a list of subsets, test for correct inclusion
(\code{is.nested()}), or strict inclusion (\code{is.strict()}).
}
\usage{
is.nested(subsets)
is.strict(subsets)
subsets.fun(n, levels = 4, prob = 0.7)
}
\arguments{
  \item{subsets}{In \code{is.nested()}, a list of subsets to be tested}
  \item{n}{Number of observations in the lowest level (ie level 1, the fastest
    code)}
  \item{levels}{Number of levels}
  \item{prob}{Probability of choosing an observation at level
  \eqn{n+1} given that there is one at the same place at level \eqn{n}}
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin (\code{subsets.fun()}); Peter Dalgaard (via R-help)}
\examples{
is.nested(subsets.fun(20))  # Should be TRUE

data(toyapps)
stopifnot(is.nested(subsets.toy))
}
\keyword{array}
