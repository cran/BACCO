\name{is.consistent}
\alias{is.consistent}
\title{Checks observational data for consistency with a subsets object}
\description{
Checks observational data for consistency with a subsets object: the
length of the vectors should match
}
\usage{
is.consistent(subsets, z)
}
\arguments{
  \item{subsets}{A subsets object}
  \item{z}{Data}    
}
\value{
Returns \code{TRUE} or \code{FALSE} depending on whether z is consistent
with subsets.
}
\references{
 M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{is.nested}}}
\examples{
data(toyapps)
stopifnot(is.consistent(subsets.toy,z.toy))

z.toy[[4]] <- 1:6
is.consistent(subsets.toy,z.toy)
}
\keyword{array}
