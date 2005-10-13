\name{latin.hypercube}
\alias{latin.hypercube}
\title{Latin hypercube design matrix}
\description{
  Gives a Latin hypercube design matrix with an arbitrary number of
  points in an arbitrary number of dimensions.  The \code{toy} dataset
  was generated using \code{latin.hypercube()}.
}
\usage{
latin.hypercube(n, d, normalize=FALSE)
}
\arguments{
  \item{n}{Number of points}
  \item{d}{Number of dimensions}
  \item{normalize}{Boolean variable with \code{TRUE} meaning to
    normalize each column so the minimum is zero and the maximum is
    one.  If it takes its default \code{FALSE}, the points represent
    midpoints of \eqn{n}{n} equispaced intervals; the points thus have a minimum of
    \eqn{0.5/n}{0.5/n} and a maximum of \eqn{1-0.5/n}{1-0.5/n}.}
}
\author{Robin K. S. Hankin}
\examples{
#10 points, 6 dimensions:
(latin.hypercube(10,6) -> x)
plot(as.data.frame(x))
}
\keyword{models}
