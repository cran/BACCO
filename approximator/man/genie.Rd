\name{genie}
\alias{genie}
\alias{Genie}
\alias{D1.genie}
\alias{hpa.genie}
\alias{z.genie}
\alias{subsets.genie}
\alias{basis.genie}
\alias{hpa.fun.genie}
\alias{hpa.genie.start}
\alias{hpa.genie.optimal}
\docType{data}
\title{Genie datasets for approximator package}
\description{
  Genie datasets that illustrate the package.
}
\usage{
  data(genie)
  D1.genie
  hpa.genie
  z.genie
  subsets.genie
  basis.genie(x)
  hpa.fun.genie(x)
  hpa.genie.start
  hpa.genie.optimal
}
\arguments{
  \item{x}{A 4-element vector (for \code{basis.genie()}); a  19-element
    vector (for \code{hpa.fun.genie()})}
  }
\format{

  The genie example is a case with three levels.
  
  The \code{D1.genie} matrix is 36 rows of code run points,
  corresponding to the observations of the level 1 code.  It has four
  columns, one per parameter.

  \code{hpa.genie} is a hyperparameter object.  

  \code{subsets.genie} is a list of three elements.  Element \eqn{i}
  corresponds to the rows of \code{D1.genie} at which level \eqn{i} has
  been observed.
  
  \code{z.genie} is a three element list.  Each element is a vector;
  element \eqn{i} corresponds to observations of level \eqn{i}.  The
  lengths will match those of \code{subsets.genie}.
 
  Function \code{basis.genie()} is a suitable basis function.
  
  Function \code{hpa.fun.genie()} creates a hyperparameter object in a
  form suitable for passing to the other functions in the library.
}
\details{

}
\references{
  M. C. Kennedy and A. O'Hagan 2000. \dQuote{Predicting the output from
    a complex computer code when fast approximations are available}
  Biometrika, 87(1): pp1-13
  }
\author{Robin K. S. Hankin}
\examples{
data(genie)
z.genie


jj <- list(trace=100,maxit=10)

hpa.genie.level1 <- opt.1(D=D1.genie, z=z.genie, basis=basis.genie, subsets=subsets.genie, hpa.start=hpa.genie.start,control=jj)
 
hpa.genie.level2 <- opt.gt.1(level=2, D=D1.genie, z=z.genie, basis=basis.genie, subsets=subsets.genie, hpa.start=hpa.genie.level1,control=jj)

hpa.genie.level3 <- opt.gt.1(level=3, D=D1.genie, z=z.genie, basis=basis.genie, subsets=subsets.genie, hpa.start=hpa.genie.level2,control=jj) 


}
\keyword{datasets}
