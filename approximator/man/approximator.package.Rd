\name{approximator-package}
\alias{approximator-package}
\alias{approximator}
\docType{package}
\title{
Bayesian approximation of computer models when fast approximations are available
}
\description{
Implements the ideas of Kennedy and O'Hagan 2000 (see references).
}
\details{
\tabular{ll}{
Package: \tab approximator\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-01-10\cr
License: \tab GPL\cr
}

This package implements the Bayesian approximation techniques discussed
in Kennedy and O'Hagan 2000.

In its simplest form, it takes input from a \dQuote{slow} code and a
\dQuote{fast} code, each run at different points in parameter space.
The approximator package then uses both sets of model runs to infer what
the top level code would produce at a given, untried point in parameter space.

}
\author{
  Robin K. S. Hankin
  
Maintainer:  <r.hankin@noc.soton.ac.uk>
}
\references{
  R. K. S. Hankin 2005. \dQuote{Introducing BACCO, an R bundle for
    Bayesian analysis of computer code output}, Journal of Statistical
  Software, 14(16)

   M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
}
\keyword{package}
\examples{
data(toyapps)
mdash.fun(x=1:3, D1=D1.toy, subsets=subsets.toy, hpa=hpa.toy, z=z.toy, basis=basis.toy)
}
