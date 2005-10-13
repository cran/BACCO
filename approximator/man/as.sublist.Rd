\name{as.sublist}
\alias{as.sublist}
\title{Converts a level one design matrix and a subsets object into a
  list of design matrices, one for each level}
\description{
Given a level one design matrix, and a subsets object, convert into a
list of design matrices, each one of which is the design matrix for its
own level
}
\usage{
as.sublist(D1, subsets)
}
\arguments{
  \item{D1}{Design matrix for level one code}
  \item{subsets}{subsets object}
}
\references{
 M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)
as.sublist(D1=D1.toy , subsets=subsets.toy)
}
\keyword{array}
