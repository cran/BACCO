\name{hpa.fun.toy}
\alias{hpa.fun.toy}
\title{Toy example of a hyperparameter object creation function}
\description{
Creates a hyperparameter object from a vector of length 19.  Intended as
a toy example to be modified for real-world cases.
}
\usage{
hpa.fun.toy(x)
}
\arguments{
  \item{x}{Vector of length 19 that specifies the correlation scales}
}
\details{
  Elements 1-4 of x specify the sigmas for each of the four levels in
  the toy example.  Elements 5-7 specify the correlation scales for
  level 1, elements 8-10 the scales for level 2, and so on.

  Internal function \code{pdm.maker()} shows how the B matrix is
  obtained from the various elements of input argument \code{x}.  Note
  how, in this simple example, the B matrices are diagonal, but
  generalizing to non-diagonal matrices should be straightforward (if
  you can guarantee that they remain positive definite).
  
}
\value{
  \item{sigmas}{The four sigmas corresponding to the four levels}
  \item{B}{The four B matrices corresponding to the four levels}
  \item{rhos}{The three (sic) matrices corresponding to levels 1-3}
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
hpa.fun.toy(1:19)
}
\keyword{array}
