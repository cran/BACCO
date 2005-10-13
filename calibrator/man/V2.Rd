\name{V2}
\alias{V2}
\title{distance between observation points}
\description{
  distance between observation points
}
\usage{
V2(x, other = NULL,  phi)
}
\arguments{
  \item{x}{Matrix whose rows are observation points}
  \item{other}{Second matrix; if \code{NULL}, use x}
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
\seealso{\code{\link{V1}}}
\examples{
data(toys)
V2(D2.toy,other=NULL, phi=phi.toy)
V2(D2.toy,x.vec,phi=phi.toy)
}
\keyword{array}
