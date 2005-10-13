\name{toy}
\alias{toy}
\docType{data}
\title{A toy dataset}
\description{
A dataframe consisting of 10 rows and 6 columns corresponding to 10 points in 
a six-dimensional space.} 
\usage{data(toy)}
\examples{
     data(toy)
     real.relation <- function(x){sum( (1:6)*x )}
     d.noisy <- apply(toy,1,real.relation)+rnorm(nrow(toy))/10

     fish <- rep(2,6)
     A <- corr.matrix(toy,scales=fish , power=1.5)
     Ainv <- solve(A)
     x.unknown <- rep(0.5 , 6)
     interpolant(x.unknown, d.noisy, toy, Ainv, scales=fish, g=FALSE)
     real.relation(x.unknown) ##should match estimate.
}
\keyword{datasets}
