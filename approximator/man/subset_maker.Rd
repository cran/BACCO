\name{subset_maker}
\alias{subset_maker}
\title{Create a simple subset object}
\description{
Given an integer vector whose \eqn{i^\mathrm{th}}{i-th} element is the
number of runs at level \eqn{i}, return a subset object in echelon
form. 
}
\usage{
subset_maker(x)
}
\arguments{
  \item{x}{A vector of integers}
}
\details{

  In this context, \code{x} being in \dQuote{echelon form} means that
  \itemize{
    \item \code{x} is consistent in the sense of passing
    \code{is.consistent()}
    \item For each \eqn{i}, \code{x[[i]] = 1:n} for some \eqn{n}.
    }  
}
\value{
A list object suitable for use as a \code{subset} object
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{is.consistent}}, \code{\link{is.nested}}, \code{\link{is.strict}}}
\examples{
subset_maker(c(10,4,3))


is.nested(subset_maker(c(4,9,6))) #should be FALSE
is.nested(subset_maker(c(9,6,4))) #should be TRUE

}
\keyword{math}
