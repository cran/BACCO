\name{optimal.scales}
\alias{optimal.scales}
\title{Use optimization techniques to find the optimal scales}
\description{
Uses optimization techniques (either Nelder-Mead or simulated annealing)
to find the optimal scales, or roughness lengths.
}
\usage{
optimal.scales(val, scales.start, d, use.like = TRUE,  give.answers = FALSE, ...)
}
\arguments{
  \item{val}{Dataframe with rows corresponding to points at which the
          function is known}
  \item{scales.start}{Initial guess for logs of the scales.  See details
  section for explanation}
  \item{d}{vector of observations, one for each row of \code{val}}
  \item{use.like}{Boolean, with default \code{TRUE} meaning to use
    likelihood for the objective function, and \code{FALSE}
    meaning to use a leave-out-one bootstrap estimator.}
  \item{give.answers}{Boolean, with default \code{FALSE} meaning to
    return just the roughness lengths and \code{TRUE} meaning to return
    extra information as returned by \code{optim()}.}
  \item{...}{Extra parameters passed to \code{optim()}.  See examples
    for usage of this argument.}
}
\details{
  Internally, this function works with the logarithms of the roughness
  lengths, because they are inherently positive.  However, note that the
  lengths themselves must be supplied to the function, not their
  logarithms. 
}
\value{
  If \code{give.answers} takes the default value of \code{FALSE}, a
  vector of roughness lengths is returned.  If \code{TRUE}, output from
  \code{optim()} is returned directly (note that element \code{par} is
  the logarithm of the desired roughness length as the optimization
  routine operates with the logs of the lengths as detailed above) 
}
\note{
  This function is slow to evaluate because it needs to
  calculate and invert \code{A} each time it is called, because the
  scales change from call to call.
}
\references{J. Oakley 2004. \dQuote{Estimating percentiles of uncertain
    computer code outputs}.  Applied Statistics, 53(1), pp89-93.

  J. Oakley 1999. \dQuote{Bayesian uncertainty analysis for complex
    computer codes}, PhD thesis, University of Sheffield.
  }
\author{Robin K. S. Hankin}
\seealso{\code{\link{interpolant}}}
\examples{

##First, define some scales:
fish <- c(1,1,1,1,1,4)

## and a sigmasquared value:
REAL.SIGMASQ <- 0.3

## and a real relation:
real.relation <- function(x){sum( (1:6)*x )}

## Now a design matrix:
val  <- latin.hypercube(100,6)

## apply the rea relation:
d <- apply(val,1,real.relation)

## and add some suitably correlated  Gaussian noise:
A <- corr.matrix(val,scales=fish)
d.noisy <-  as.vector(rmvnorm(n=1,mean=apply(val,1,real.relation),REAL.SIGMASQ*A))

##  Now see if we can estimate the roughness lengths well.  Remember that
##  the true values are those held in vector "fish":
optimal.scales(val=val, scales.start=rep(1,6), d=d.noisy, method="SANN",control=list(trace=1000,maxit=3), give=FALSE)
}
\keyword{models}
