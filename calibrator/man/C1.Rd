\name{C1}
\alias{C1}
\title{Matrix of distances from D1 to D2}
\description{
Returns a matrix of distances from the code run points to the
augmented observation points.  A wrapper for \code{c1.fun()}.
}
\usage{
C1(D1, D2,  theta, phi)
}
\arguments{
  \item{D1}{D1}
  \item{D2}{D2}
  \item{theta}{Parameters}
  \item{phi}{Hyperparameters}
}
\references{\dQuote{Bayesian calibration of computer models}, M. C. Kennedy and
A. O'Hagan 2001.  Journal of the Royal Statistical Society B, 63(3)
pp425-464

\dQuote{Supplementary details on Bayesian calibration of computer models},
M. C. Kennedy and A. O'Hagan 2001.   Internal report, University of
Sheffield}
\author{Robin K. S. Hankin}
\seealso{\code{\link{t.fun}}}
\examples{
data(toys)
C1(D1=D1.toy, D2=D2.toy, theta=theta.toy, phi=phi.toy)
}
\keyword{array}
