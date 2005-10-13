\name{hdash.fun}
\alias{hdash.fun}
\title{Hdash}
\description{
Returns the thing at the top of page 6
}
\usage{
hdash.fun(x, hpa, basis)
}

\arguments{
  \item{x}{Point in question}
  \item{hpa}{Hyperparameter object}
  \item{basis}{Basis functions}
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)
hdash.fun(x=1:3 , hpa=hpa.toy,basis=basis.toy)

uu <- rbind(1:3,1,1:3,3:1)
rownames(uu) <- paste("uu",1:4,sep="_")
hdash.fun(x=uu, hpa=hpa.toy,basis=basis.toy)

}
\keyword{array}
