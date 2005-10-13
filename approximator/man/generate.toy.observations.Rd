\name{generate.toy.observations}
\alias{generate.toy.observations}
\title{Er, generate toy observations}
\description{
Generates toy observations on four levels using either internal
(unknown) parameters and hyperparameters, or user-supplied versions.
}
\usage{
generate.toy.observations(D1, subsets, basis.fun, hpa = NULL, betas = NULL,
export.truth = FALSE)
}
\arguments{
  \item{D1}{Design matrix for level 1 code}
  \item{subsets}{Subset object}
  \item{basis.fun}{Basis function}
  \item{hpa}{Hyperparameter object.  If \code{NULL}, use the  internal
    (true but unknown) hyperparameter object}
  \item{betas}{Regression coefficients.  If \code{NULL}, use the  internal
    (true but unknown) regression coefficients}
  \item{export.truth}{Boolean, with default \code{FALSE} meaning to
    return synthetic observations and \code{TRUE} meaning to return the
    actual hyperparameters and coefficients.}
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)
generate.toy.observations(D1=D1.toy, subsets=subsets.toy, basis.fun=basis.toy)
}
\keyword{array}
