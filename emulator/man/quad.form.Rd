\name{quad.form}
\alias{quad.form}
\alias{quad.form.inv}
\alias{quad.3form}
\title{Evaluate a quadratic form efficiently}
\description{
  Given a square matrix \code{M} of size \eqn{n\times n}{n*n}, and a matrix
  \eqn{x} of size \eqn{n\times p}{n*p} (or a vector of length
  \eqn{n}), evaluate \eqn{x^TMx}{t(x) \%*\% M \%*\% x} in an
  efficient manner.

  Function \code{quad.form.inv()} returns \eqn{x^TM^{-1}x}{t(x) \%*\%
    solve(M) \%*\% x} using an efficient method that avoids inverting
    \eqn{M}.

  Function \code{quad.3form()} returns \eqn{l^TMr}{t(l) \%*'\% M \%*\%
  r} using nested calls to \code{crossprod()}.  It's no faster than
  calling \code{crossprod()} directly, but makes code neater (IMHO)
}
\usage{
quad.form(M, x, chol=FALSE)
quad.form.inv(M, x)
quad.3form(M,left,right)
}
\arguments{
  \item{M}{Square matrix of size \eqn{n\times n}{n*n}}
  \item{x}{Matrix of size \eqn{n\times p}{n*p}, or vector of length \eqn{n}}
  \item{chol}{Boolean, with \code{TRUE} meaning to interpret
    argument \code{M} as the lower triangular Cholesky decomposition
    of the quadratic form.  Remember that 
    \code{M.lower \%*\% M.upper == M},  and \code{chol()} returns the
    upper triangular matrix, so one needs to use the transpose
    \code{t(chol(M))}. }
  \item{left,right}{In function \code{quad.3form()}, matrices with
    \eqn{n} rows and arbitrary number of columns}
}
\details{
  The \dQuote{meat} of \code{quad.form()} for \code{chol=FALSE} is just
  \code{crossprod(crossprod(M, x), x)}, and that of
  \code{quad.form.inv()} is \code{crossprod(x, solve(M, x))}.
  
  If the Cholesky decomposition of \code{M} is available, then calling
  with \code{chol=TRUE} and supplying \code{M.upper} should generally be
  faster (for large matrices) than calling with \code{chol=FALSE} and
  using \code{M} directly.  The time saving is negligible for matrices
  smaller than about \eqn{50\times 50}{50*50}, even if the overhead of
  computing \code{M.upper} is ignored.

}
\author{Robin K. S. Hankin}
\note{These functions are used extensively in the emulator and
  calibrator packages' R code in the interests
  of speed.  They are not really intended for the end user.
}
\seealso{\code{\link{optimize}}}
\examples{
jj <- matrix(rnorm(80),20,4)
M <- crossprod(jj,jj)
M.lower <- t(chol(M))
x <- matrix(rnorm(8),4,2)

jj.1 <- t(x) \%*\% M \%*\% x
jj.2 <- quad.form(M,x)
jj.3 <- quad.form(M.lower,x,chol=TRUE)

print(jj.1)
print(jj.2)
print(jj.3)

## Now consider accuracy:
quad.form(solve(M),x) - quad.form.inv(M,x)  # should be zero
}
\keyword{array}
