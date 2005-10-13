\name{create.new.toy.datasets}
\alias{create.new.toy.datasets}
\docType{data}
\title{Create new toy datasets}
\description{
  Creates new toy datasets, by sampling from an explicitly specified
  multivariate Gaussian distribution whose covariance matrix is that
  required for a Gaussian process.
}

\usage{create.new.toy.datasets(N=8,n=5,num.of.var.params=2,
num.of.code.params=3, set.seed.to.zero=FALSE,normalize=FALSE)}
\arguments{
  \item{N}{Number of code runs}
  \item{n}{number of observation points}
  \item{num.of.var.params}{number of parameters}
  \item{num.of.code.params}{number of code parameters}
  \item{set.seed.to.zero}{Boolean, with \code{TRUE} meaning to set the
    seed to zero}
  \item{normalize}{Boolean, with \code{TRUE} meaning to return columns
    whose min is zero and max is one.  See \code{latin.hypercube} in
    \code{library(emulator)} for details.}
  }
\value{
    Returns a list of seven elements:
    \item{x.star.toy}{}
    \item{t.vec.toy}{}
    \item{y.toy}{}
    \item{z.toy}{}
    \item{d.toy}{}
    \item{D1.toy}{}
    \item{D2.toy}{}

    these being the toy datasets that are needed for the examples.
    Datasets created with \code{create.new.toy.datasets()} will work
    with the examples whatever the value of \code{n} and \code{N}
    (although \code{num.of.var.params} and \code{num.of.code.params}
    must take their default values).    
}
\note{
  Because function \code{create.new.toy.datasets()} calls
    \code{computer.model()} and \code{reality()}, the datasets returned
    are drawn from a multivariate Gaussian distribution which \strong{is} a
    Gaussian process.
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
\seealso{
  \code{\link{toys}}, \code{\link{reality}}, \code{\link{latin.hypercube}}
}
\examples{
create.new.toy.datasets()
create.new.toy.datasets(N=40)
}
\keyword{datasets}
