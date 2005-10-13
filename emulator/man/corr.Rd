\name{corr}
\alias{corr}
\alias{corr.matrix}
\title{correlation function for calculating A}
\description{
calculates the correlation function between two points in parameter
space, thus determining the correlation matrix A.
}
\usage{
corr(x1, x2, scales=NULL , pos.def.matrix=NULL,
power=2,coords="cartesian", spherical.distance.function=NULL)
corr.matrix(xold, yold, use.neatversion=TRUE, distance.function=corr, ...)
}
\arguments{
  \item{x1}{First point}
  \item{x2}{Second point}
  \item{scales}{Vector specifying the diagonal elements of B (see below)}
  \item{pos.def.matrix}{Positive definite matrix to be used by
    \code{corr.matrix()} for B.
    Exactly one of \code{scales} and \code{pos.definite.matrix} should
    be specified.  Supplying \code{scales} specifies the diagonal
    elements of B (off diagonal elements are set to zero); supply
     \code{pos.definite.matrix} in the general case.  Note that
    neither \code{corr()} nor \code{corr.matrix()} test for positive
    definiteness.}
  \item{power}{The power to use in the exponent of the exponential.  In
    the case of a diagonal \code{pos.def.matrix} (or,
    equivalently, using \code{scales}), the default value of 2 gives the
    standard metric, viz \eqn{e^{-(x-x')^2}}{exp(-(x-x')^2)}.

    For the general case of nondiagonal \code{pos.def.matrix},
    the metric is \eqn{e^{-\left| (x-x')^TB(x-x')\right|^{\rm
    power/2}}}{exp(-abs( (x-x')^T B (x-x') )^(power/2))}, thus reducing
    to the standard case for \code{power}=2.}
  \item{coords}{In function \code{corr()}, a character string, with
    default \dQuote{cartesian} meaning to interpret the elements of
    \code{x1} (and \code{x2}) as coordinates in Cartesian space.  The
    only other acceptable value is currently \dQuote{spherical}, which
    means to interpret the first element of \code{x1} as row number, and
    the second element as column number, on a spherical computational
    grid (such as used  by climate model Goldstein; see package
    \code{goldstein} for an example of this option in use).
    Spherical geometry is then used to calculate the geotetic (great
    circle) distance between point \code{x1} and \code{x2}, with
    function \code{gcd()}.}
  \item{distance.function}{Function to be used to calculate distances.
    Defaults to \code{corr()}.}
  \item{xold}{data frame each row of which is an evaluated point}
  \item{yold}{(optional) data frame each row of which is an evaluated
    point.  If missing, use \code{xold}}
  \item{use.neatversion}{Boolean variable passed to \code{corr.matrix()} with
    \code{TRUE} meaning to use the iterated \code{apply()} method
    and \code{FALSE} meaning to use two nested loops.}
  \item{spherical.distance.function}{In \code{corr}, a function to
    determine the distance between two points; used if
    \code{coords}=\dQuote{spherical}. A good one to choose is
    \code{gcd()} (that is, Great Circle Distance) of the goldstein library.} 
  \item{...}{extra arguments given to \code{corr.matrix()} that are
    passed on to \code{corr()} (such as \code{scales} and \code{power}
    if \code{distance.function=corr}).}
}

\details{
  Evaluates a slight generalization of Oakley's equation 2.12 for the
  correlation between \eqn{\eta(x)}{eta(x)} and \eqn{\eta(x')}{eta(x')}:
  \eqn{e^{-\left| (x-x')^TB(x-x')\right|^{\rm power/2}}}{exp(|(x-x')^T B
    (x-x')|^({power/2))}}.
  This reduces to Oakley's form if \code{power=2}.
  }
\value{
Returns the correlation function
}
\references{
  J. Oakley 1999. \dQuote{Bayesian uncertainty analysis for complex
    computer codes}, PhD thesis, University of Sheffield.

  J. Oakley and A. O'Hagan, 2002. \dQuote{Bayesian Inference for the
  Uncertainty Distribution of Computer Model Outputs}, Biometrika
  89(4), pp769-784

  R. K. S. Hankin 2005. \dQuote{Introducing BACCO, an R bundle for
    Bayesian analysis of computer code output}, Journal of Statistical
  Software, 14(16)
}
\author{Robin K. S. Hankin}
\examples{
corr(1:10,10:1,scales=rep(1,10), power=2)
corr(1:10,10:1,pos.def.matrix=0.1+diag(10),power=2)
x <- latin.hypercube(4,7)  #4 points in 7-dimensional space
corr.matrix(x,scales=rep(1,7),power=1.5)

x[1,1] <- 100
corr.matrix(x,scales=rep(1,7), power=1.5)


# note that all the first row and first column apart from the [1,1]th
# element is zero (or very nearly so) because point number 1 is now very
# far from the other points.

#to use just a single dimension, remember to use the drop=FALSE argument:
corr.matrix(x[,1,drop=FALSE],scales=rep(1,1),power=1.5)


# For problems in 1D, coerce the independent variable to a matrix:
m <- c(0.2 , 0.4,0.403,0.9)
corr.matrix(cbind(m),scales=1)


# now use a non-default value for distance.function.
# Function f() below taken from Oakley's thesis page 12,
# equation 2.10:

f <- function(x,y,theta){
  d <- sum(abs(x-y))
  if(d >= theta){
    return(0)
  }else{
    return(1-d/theta)
  }
}

corr.matrix(xold=x, distance.function=f , theta=4)
 # Note the first row and first column is a single 1 and 3 zeroes
 # (because the first point, viz x[1,], is "far" from the other ponts).



}

\keyword{models}
