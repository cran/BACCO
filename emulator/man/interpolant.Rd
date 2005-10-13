\name{interpolant}
\alias{interpolant}
\alias{interpolant.quick}
\title{Interpolates between known points using Bayesian estimation}
\description{Calculates the a postiori distribution of results at a
  point using the techniques outlined by Oakley.  This function is the
  primary function of the package.  Function \code{interpolant.quick()}
  gives the expectation of the emulator at a set of points, and function
  \code{interpolant()} gives the expectation and other information (such
  as the variance) at a single point.  
}
\usage{
interpolant(x, d, xold, Ainv=NULL, A=NULL, use.Ainv=TRUE, scales=NULL, pos.def.matrix=NULL,
func=regressor.basis,  power=2, give.full.list = FALSE)
interpolant.quick(x, d, xold, Ainv, scales=NULL,
pos.def.matrix=NULL, func=regressor.basis, give.Z = FALSE, power=2)
}
\arguments{
  \item{x}{Point(s)  at which estimation is desired.  For
    \code{interpolant.quick()}, argument \code{x} is a data frame
    and an expectation is given for each row.}
  \item{d}{vector of observations, one for each row of \code{xold}}
  \item{xold}{Data frame with rows corresponding to points at which the
    function is known}
  \item{A}{Correlation matrix \code{A}.  If not given, it is calculated.}
  \item{Ainv}{Inverse of correlation matrix \code{A}.  Required by
    \code{interpolant.quick()}.  In \code{interpolant()},
    using the default value of \code{NULL} results in \code{Ainv} being
    calculated explicitly (which may be slow: see next argument for more
    details).}
  \item{use.Ainv}{Boolean, with default \code{TRUE} meaning to use the
    inverse matrix \code{Ainv} (and, if necessary, calculate it using
    \code{solve(.)}).  This requires the not inconsiderable overhead of
    inverting a matrix.   If, however, \code{Ainv} is available, using
    the default option is \emph{much} faster than setting
    \code{use.Ainv=FALSE}; see below.
    
    If \code{FALSE}, function \code{interpolant()} does not use
    \code{Ainv}, but makes extensive use of \code{solve(A,x)} (mostly in
    the form of \code{quad.form.inv()} calls).  This option avoids the
    overhead of inverting a matrix, but has non-negligible marginal costs.

    If \code{Ainv} is not available, there is little to choose, in terms
    of execution time, between calculating it explicitly (that is,
    setting \code{use.Ainv=TRUE}) and using \code{solve(A,x)} (ie
    \code{use.Ainv=TRUE}).
    
    \strong{Note:} if \code{Ainv} is given to the function, but
    \code{use.Ainv} is \code{FALSE}, the code will do as requested and use
    the slow \code{solve(A,x)}, which is probably not what you want.  }
    
  \item{func}{Function used to determine basis vectors, defaulting
    to \code{regressor.basis}} if not given.
  \item{give.full.list}{In \code{interpolant()}, Boolean variable with
    \code{TRUE} meaning to return the whole list of a postiori
    parameters as detailed on pp12-15 of Oakley, and default \code{FALSE}
    meaning to return just the best estimate.}
  \item{scales}{Vector of \dQuote{roughness} lengths used to calculate
    \code{t(x)}.  Note that \code{scales} is needed twice: once to calculate
    \code{Ainv} and once to calculate \code{t(x)} inside
    \code{interpolant} (which is determined by calling \code{corr}
    inside an \code{apply()} loop).  A good place to start might be
    \code{scales=rep(1,ncol(xold))}.}
  \item{pos.def.matrix}{A positive definite matrix that is used if
    \code{scales} is not supplied.  Note that precisely one of
    \code{scales} and \code{pos.def.matrix} must be supplied.}
  \item{give.Z}{In function \code{interpolant.quick()}, Boolean variable
    with \code{TRUE} meaning to return the best estimate and the error,
    and default \code{FALSE} meaning to return just the best estimate.}
  \item{power}{The power to use in the correlation function.  Default of
    2 corresponds to Oakley's equation 2.12 (see \code{?corr} for more
    details).   Note that if \code{Ainv} is supplied, then the
    \code{power} used to calculate \code{Ainv} must match the power used
    in \code{interpolant()}
  }
  }
\value{
  If \code{give.full.list} is TRUE, a list is return with components
  \item{betahat}{Standard MLE of the (linear) fit, given the
    observations}
  \item{prior}{Estimate for the prior}
  \item{sigmahat.square}{A postiori estimate for variance}
  \item{mstar.star}{A postiori expectation}
  \item{cstar}{a priori correlation of a point with itself}
    \item{cstar.star}{A postiori correlation of a point with itself}
  \item{Z}{Standard deviation (although the distribution is actually a
    t-distribution with n-q degrees of freedom)}
}
\references{
  J. Oakley 2004. \dQuote{Estimating percentiles of uncertain
    computer code outputs}.  Applied Statistics, 53(1), pp89-93.

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
\seealso{\code{\link{makeinputfiles}}}
\examples{
# example has 10 observations on 6 dimensions.
# function is just sum( (1:6)*x) where x=c(x_1, ... , x_2)

data(toy)
val <- toy
real.relation <- function(x){sum( (0:6)*x )}
H <- regressor.multi(val)
d <- apply(H,1,real.relation)

fish <- rep(1,6)
fish[6] <- 4

A <- corr.matrix(val,scales=fish, power=2)
Ainv <- solve(A)

# now add some suitably correlated noise to d:
d.noisy <-  as.vector(rmvnorm(n=1, mean=d, 0.1*A))
names(d.noisy) <- names(d)

# First try a value at which we know the answer (the first row of val):
x.known <- as.vector(val[1,])
bayes.known <- interpolant(x.known, d, val, Ainv=Ainv, scales=fish, g=FALSE)
print("error:")
print(d[1]-bayes.known)


# Now try the same value, but with noisy data:
print("error:")
print(d.noisy[1]-interpolant(x.known, d.noisy, val, Ainv=Ainv, scales=fish, g=FALSE))

#And now one we don't know:
x.unknown <- rep(0.5 , 6)
bayes.unknown <- interpolant(x.unknown, d.noisy, val, scales=fish, Ainv=Ainv,g=TRUE)

## [   compare with the "true" value of sum(0.5*0:6) = 10.5   ]


## (To find the best correlation lengths, use optimal.scales())



 # Now we use the SAME dataset but a different set of basis functions.
 # Here, we use the functional dependence of
 # "A+B*(x[1]>0.5)+C*(x[2]>0.5)+...+F*(x[6]>0.5)".
 # Thus the basis functions will be c(1,x>0.5).
 # The coefficients will again be 1:6.

       # Basis functions:
f <- function(x){c(1,x>0.5)}
       # (other examples might be
       # something like  "f <- function(x){c(1,x>0.5,x[1]^2)}"

       # now create the data
real.relation2 <- function(x){sum( (0:6)*f(x) )}
d2 <- apply(val,1,real.relation2)

       # Define a point at which the function's behaviour is not known:
x.unknown2 <- rep(1,6)
       # Thus real.relation2(x.unknown2) is sum(1:6)=21

       # Now try the emulator:
interpolant(x.unknown2, d2, val, Ainv=Ainv, scales=fish, g=TRUE)$mstar.star
       # Heh, it got it wrong!  (we know that it should be 21)


       # Now try it with the correct basis functions:
interpolant(x.unknown2, d2, val, Ainv=Ainv,scales=fish, func=f,g=TRUE)$mstar.star
       # That's more like it.

       # We can tell that the coefficients are right by:
betahat.fun(val,Ainv,d2,func=f)
       # Giving c(0:6), as expected.

       # It's interesting to note that using the *wrong* basis functions
       # gives the *correct* answer when evaluated at a known point:
interpolant(val[1,], d2, val, Ainv=Ainv,scales=fish, g=TRUE)$mstar.star
real.relation2(val[1,])
       # Which should agree.


       # Now look at Z.  Define a function Z() which determines the
       # standard deviation at a point near a known point.
Z <- function(o) {
    x <- x.known 
    x[1] <- x[1]+ o
    interpolant(x, d.noisy, val, Ainv=Ainv, scales=fish, g=TRUE)$Z.oak
  } 

Z(0)       #should be zero because we know the answer (this is just Z at x.known)
Z(0.1)     #nonzero error.


  ## interpolant.quick() should  give the same results faster, but one
  ##   needs a matrix:
u <- rbind(x.known,x.unknown)
interpolant.quick(u, d.noisy, val, scales=fish, Ainv=Ainv,g=TRUE)

data(results.table)
data(expert.estimates)


       # Decide which column we are interested in:
output.col <- 26

       #
wanted.cols <- c(2:9,12:19)

       # Decide how many to keep;
       # 30-40 is about the most we can handle:
wanted.row <- 1:27

       # Values to use are the ones that appear in goin.test2.comments:
val <- results.table[wanted.row , wanted.cols]

       # Now normalize val so that 0<results.table[,i]<1 for all i:

normalize <- function(x){(x-mins)/(maxes-mins)}
unnormalize <- function(x){mins + (maxes-mins)*x}

mins  <- expert.estimates$low 
maxes <- expert.estimates$high
jj <- t(apply(val,1,normalize))

jj <- as.data.frame(jj) 
names(jj) <- names(val)
val <- jj


       ## Answer is the 19th (or 20th or ... or 26th)
d  <- results.table[wanted.row ,  output.col]

A <- corr.matrix(val,scales=rep(1,ncol(val)), power=1.5)
Ainv <-  solve(A)


scales.optim <- c( -2.917, -4.954, -3.354, 2.377, -2.457, -1.934, -3.395,
-0.444, -1.448, -3.075, -0.052, -2.890, -2.832, -2.322, -3.092, -1.786)

print("and plot points used in optimization:")
d.observed <- results.table[ , output.col]

A <- corr.matrix(val,scales=scales.optim, power=1.5)
Ainv <- solve(A)

print("now plot all points:")
design.normalized <- as.matrix(t(apply(results.table[,wanted.cols],1,normalize)))
d.predicted <- interpolant.quick(design.normalized , d , val , Ainv=Ainv,
scales=scales.optim, power=1.5)
jj <- range(c(d.observed,d.predicted))
par(pty="s")
plot(d.observed, d.predicted, pch=16, asp=1,
xlim=jj,ylim=jj,
xlab=expression(paste(temperature," (",{}^o,C,"), model"   )),
ylab=expression(paste(temperature," (",{}^o,C,"), emulator"))
)
abline(0,1)
}
\keyword{models}
