\name{cov.p5.supp}
\alias{cov.p5.supp}
\title{Covariance function for posterior distribution of z}
\description{
Covariance function for posterior distribution of \eqn{z(\cdot)}{z(.)}
conditional on estimated hyperparameters and calibration parameters
\eqn{\theta}{theta}. 
}
\usage{
cov.p5.supp(x, xdash, theta, d, D1, D2, H1, H2, phi)
}
\arguments{
  \item{x}{first point, or a matrix whose rows are the points of interest}
  \item{xdash}{second point, or a matrix whose rows are the points of interest}
  \item{theta}{Parameters}
  \item{d}{Observed values}
  \item{D1}{Code run design matrix}
  \item{D2}{Observation points of real process}
  \item{H1}{Basis function for \code{D1}}
  \item{H2}{Basis function for \code{D2}}
  \item{phi}{Hyperparameters}
}
\details{
Evaluates the covariance function: the last formula on page 5 of the supplement
}
\value{
Returns a matrix of covariances
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
\note{
  May return the transpose of the desired object}
  }
  \examples{
data(toys)
x <- rbind(x.toy,x.toy+1,x.toy,x.toy,x.toy)
rownames(x) <- letters[1:5]
xdash <- rbind(x*2,x.toy)
rownames(xdash) <- LETTERS[1:6]

cov.p5.supp(x=x,xdash=xdash,theta=theta.toy,d=d.toy,D1=D1.toy,D2=D2.toy,H1=H1.toy,H2=H2.toy,phi=phi.toy)


phi.true <- phi.true.toy(phi=phi.toy)
cov.p5.supp(x=x,xdash=xdash,theta=theta.toy,d=d.toy,D1=D1.toy,D2=D2.toy,H1=H1.toy,H2=H2.toy,phi=phi.true)
}
\keyword{array}
