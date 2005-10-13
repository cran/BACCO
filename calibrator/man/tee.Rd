\name{tee}
\alias{tee}
\title{tee: distances to a point}
\description{
  Returns a vector whose elements are the \dQuote{distances} from a point
  to the observations and code run points.
}
\usage{
tee(x, theta, D1, D2, phi)
}
\arguments{
  \item{x}{Point from which distances are calculated}
  \item{theta}{Value of parameters}
  \item{D1}{Matrix of code run points}
  \item{D2}{Matrix of observation points}
  \item{phi}{Hyperparameters}
}
\details{
  This function returns the column vector just below equation 10 of
  KOH2001 (NB: not equation 10 of the supplement!)
}
\references{
  M. C. Kennedy and A. O'Hagan 2001. \dQuote{Bayesian
    calibration of computer models}.  Journal of the Royal Statistical
  Society B, 63(3) pp425-464
  
  M. C. Kennedy and A. O'Hagan 2001.  \dQuote{Supplementary details on
    Bayesian calibration of computer models}, Internal report, University
  of Sheffield.  Available at \url{http://www.shef.ac.uk/~st1ao/ps/calsup.ps}
  
  R. K. S. Hankin 2005. \dQuote{Introducing BACCO, an R bundle for
    Bayesian analysis of computer code output}, Journal of Statistical
  Software, 14(16)
}
\author{Robin K. S. Hankin}
\seealso{
  \code{\link{Ez.eqn9.supp}}
}
\examples{
data(toys)
tee(x=x.toy, theta=theta.toy, D1=D1.toy, D2=D2.toy,
phi=phi.toy)
}
\keyword{array}
