\name{etahat}
\alias{etahat}
\title{Expectation of computer output}
\description{
Returns the apostiori expectation of the computer program at a particular point
with a particular set of parameters, given the code output.
}
\usage{
etahat(D1, D2, H1, y, E.theta, extractor, phi)
}

\arguments{
  \item{D1}{Matrix of code observation points and parameters}
  \item{D2}{Matrix of field observation points}
  \item{H1}{Basis functions}
  \item{y}{Point at which code's expectation is calculated (for
    \code{etahat.individual})}
  \item{E.theta}{expectation wrt theta}
  \item{extractor}{Extractor function}
  \item{theta}{Parameters}
  \item{phi}{Hyperparameters}
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
\seealso{\code{\link{p.page4}}}
\examples{
data(toys)

etahat(D1=D1.toy, D2=D2.toy, H1=H1.toy, y=y.toy, E.theta=E.theta.toy,
extractor=extractor.toy, phi=phi.toy)
}
\keyword{array}
