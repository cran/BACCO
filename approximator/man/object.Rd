\name{object}
\alias{object}
\alias{opt.1}
\alias{opt.gt.1}
\title{Optimization of posterior likelihood of hyperparameters}
\description{
  Returns the likelihood of a set of hyperparameters given the data.
  Functions \code{opt1()} and \code{opt.gt.1()} find hyperparameters
  that maximize the relevant likelihood for level 1 and higher levels respectively.
}
\usage{
object(level, D, z, basis, subsets, hpa)
opt.1(D, z, basis, subsets, hpa.start, give.answers=FALSE, ...)
opt.gt.1(level, D, z, basis, subsets, hpa.start, give.answers=FALSE, ...)
}
\arguments{
  \item{level}{level}
  \item{D}{Design matrix for top-level code}
  \item{z}{Data}
  \item{basis}{Basis function}
  \item{subsets}{subsets object}
  \item{hpa}{hyperparameter object}
  \item{hpa.start}{Starting value for hyperparameter object}
  \item{give.answers}{Boolean, with default \code{FALSE} meaning to
    return just the point estimate, and \code{TRUE} meaning to return
    extra information from the call  to \code{optim()}.}
  \item{...}{Extra arguments passed to \code{optim()}}
}
\details{
This function is the object function used in toy optimizers
\code{optimal.hpa()}. 
}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(toyapps)
object(level=4, D=D1.toy , z=z.toy,basis=basis.toy,
   subsets=subsets.toy, hpa=hpa.fun.toy(1:19))
object(level=4, D=D1.toy , z=z.toy,basis=basis.toy,
   subsets=subsets.toy, hpa=hpa.fun.toy(3+(1:19)))
} 
\keyword{array}