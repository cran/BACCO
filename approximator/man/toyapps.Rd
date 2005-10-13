\name{toyapps}
\alias{D1.toy}
\alias{hpa.toy}
\alias{z.toy}
\alias{subsets.toy}
\alias{betas.toy}
\docType{data}
\title{Toy datasets for approximator package}
\description{
  Toy datasets that illustrate the package.
}
\usage{
  data(toyapps)
  D1.toy
  hpa.toy
  z.toy
  subsets.toy
  betas.toy
}
\format{

  The toy example is a case with four levels.
  
  The \code{D1.toy} matrix is 20 rows of code run points, corresponding
  to the observations of the level 1 code.  It has three columns, one
  per parameter.

  \code{hpa.toy} is a hyperparameter object.  It is a list of three
  elements: \code{sigmas}, \code{B}, and \code{rhos}.

  \code{subsets.toy} is a list of four elements.  Element \eqn{i}
  corresponds to the rows of \code{D1.toy} at which level \eqn{i} has
  been observed.
  
  \code{z.toy} is a four element list.  Each element is a vector;
  element \eqn{i} corresponds to obsevations of level \eqn{i}.  The
  lengths will match those of \code{subsets.toy}.

  \code{betas.toy} is a matrix of coefficients.
  
  \strong{Brief description of toy functions fully documented under their own manpage}
  
  Function \code{generate.toy.observations()} creates new toy datasets
  with any number of observations and code runs.

  Function \code{basis.toy()} is an example of a basis function
  
  Function \code{hpa.fun.toy()} creates a hyperparameter object such as
  \code{phi.toy} in a form suitable for passing to the other functions
  in the library.

  \strong{See the helpfiles listed in the \dQuote{see also} section
  below}
  
}
\details{
  All toy datasets are documented here.  There are also several toy
  functions that are needed for a toy problem; these are documented
  separately (they are too diverse to document fully in a single
  manpage).  Nevertheless a terse summary  for each toy function
  is provided on this page.  All toy functions in the package are listed
  under \dQuote{See Also}.
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)

is.consistent(subsets.toy , z.toy)

generate.toy.observations(D1.toy, subsets.toy, basis.toy, hpa.toy, betas.toy)

}
\keyword{datasets}
