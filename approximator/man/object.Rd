\name{object}
\alias{object}
\alias{opt.1}
\alias{opt.gt.1}
\title{Optimization of posterior likelihood of hyperparameters}
\description{
  Returns the likelihood of a set of hyperparameters given the data.
  Functions \code{opt1()} and \code{opt.gt.1()} find hyperparameters
  that maximize the relevant likelihood for level 1 and higher levels
  respectively.   Function \code{object()} returns the expression given
  by equation 9 in KOH2000, which is minimized \code{opt1()} and
  \code{opt.gt.1()}.
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
  \item{...}{Extra arguments passed to \code{optim()}.  A common one
    would be \code{control=list(trace=100)}.}
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
\seealso{\code{\link{genie}}}  
\examples{
data(toyapps)
object(level=4, D=D1.toy , z=z.toy,basis=basis.toy,
   subsets=subsets.toy, hpa=hpa.fun.toy(1:19))
object(level=4, D=D1.toy , z=z.toy,basis=basis.toy,
   subsets=subsets.toy, hpa=hpa.fun.toy(3+(1:19)))


# Now a little example of finding optimal hyperpameters in the toy case
# (a bigger example is given on the genie help page)
jj <- list(trace=100,maxit=10)

hpa.toy.level1 <- opt.1(D=D1.toy, z=z.toy, basis=basis.toy, subsets=subsets.toy, hpa.start=hpa.toy,control=jj)

hpa.toy.level2 <- opt.gt.1(level=2, D=D1.toy, z=z.toy, basis=basis.toy, subsets=subsets.toy, hpa.start=hpa.toy.level1,control=jj)

hpa.toy.level3 <- opt.gt.1(level=3, D=D1.toy, z=z.toy, basis=basis.toy, subsets=subsets.toy, hpa.start=hpa.toy.level2,control=jj) 

hpa.toy.level4 <- opt.gt.1(level=4, D=D1.toy, z=z.toy, basis=basis.toy, subsets=subsets.toy, hpa.start=hpa.toy.level3,control=jj)

} 
\keyword{array}
