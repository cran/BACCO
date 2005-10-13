\name{Pi}
\alias{Pi}
\title{Kennedy's Pi notation}
\description{
Evaluates Kennedy's \eqn{\prod}{Pi} product
}
\usage{
Pi(hpa, i, j)
}
\arguments{
  \item{hpa}{Hyperparameter object}
  \item{i}{subscript}
  \item{j}{superscript}
}
\details{
This function evaluates Kennedy's \eqn{\prod}{Pi} product, but with the
additional feature that \eqn{\prod_i^j=0}{Pi_i^j=0} if \eqn{i>j+1}.
This seems to work in practice.
}
\references{
 M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
}
\author{Robin K. S. Hankin}
\examples{
data(toyapps)
Pi(hpa.toy,1,2)
Pi(hpa.toy,2,2)
Pi(hpa.toy,3,2)
Pi(hpa.toy,4,2)
}
\keyword{array}
