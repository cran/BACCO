\name{A}
\alias{A}
\title{Matrix of correlations between two sets of points}
\description{
  Returns the matrix of correlations of code output at each level
  evaluated at points on the design matrix.
}
\usage{
A(level, Di, Dj, hpa)
}
\arguments{
  \item{level}{The level.  This enters via the correlation scales}
  \item{Di}{First set of points}
  \item{Dj}{Second set of points}
  \item{hpa}{Hyperparameter object}
}
\details{
  This is essentially a convenient wrapper for function
  \code{corr.matrix}.   It is not really intended for the end user.  
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
}
\author{Robin K. S. Hankin}
\seealso{\code{\link[emulator]{corr.matrix}},\code{\link{c_fun}}}
\examples{
data(toyapps)
D2 <-  D1.toy[subsets.toy[[2]],]
D3 <-  D1.toy[subsets.toy[[3]],]

A(1,D2,D3,hpa.toy)
A(2,D2,D3,hpa.toy)

}
\keyword{array}
