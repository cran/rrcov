\name{rrcov.control}
\alias{rrcov.control}
\alias{rrcov.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control object for the estimation parameters }
\description{
     Useful for passing the estimation options as parameters to the estimation functions
}
\usage{
rrcov.control(alpha=1/2, nsamp=500, seed=NULL, tolSolve=10e-14,
    trace=FALSE, use.correction=TRUE, adjust=FALSE,
    r = 0.45, arp = 0.05, eps=1e-3, maxiter=120)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{ This parameter controls the size of the subsets over 
   which the determinant is minimized, i.e. \code{alpha*n} 
   observations are used for computing the determinant. Allowed 
   values are between 0.5 and 1 and the default is 0.5. }
  \item{nsamp}{number of subsets used for initial estimates or \code{"best"} 
    or \code{"exact"}. Default is \code{nsamp = 500}. 
    If \code{nsamp="best"} exhaustive enumeration is done, as far as 
    the number of trials do not exceed 5000. If \code{nsamp="exact"} 
    exhaustive enumeration will be attempted however many samples 
    are needed. In this case a warning message will be displayed 
    saying that the computation can take a very long time. }
  \item{seed}{ starting value for random generator. Default is \code{seed = NULL}}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.}
  \item{trace}{ whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{use.correction}{ whether to use finite sample correction factors. Default is \code{use.correction=TRUE}}
  \item{adjust}{ whether to perform intercept adjustment at each step. This could be quite 
time consuming, therefore the default is \code{adjust = FALSE} }
  \item{r}{ M-estimates: breakdown point, i.e. the fraction of contaminated data. The default is 0.45}
  \item{arp}{ M-estimates: asypmthotic rejection point, i.e. the fraction of points 
    receiving zero weights. The default is 0.001}
  \item{eps}{ M-estimates: the relaive precision of the solution. The default is 1e-3}
  \item{maxiter}{ M-estimates: maximum number of iterations for the computation of the M-estimates. The default is 120}
}
\details{
For details about the estimation options see the corresponding estimation functions.
}
\value{
  A list with components, as the parameters passed  by the invocation

}
\references{ 
}

% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{}
\examples{
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])

ctrl <- rrcov.control(alpha=0.75, trace=TRUE)
covMcd(hbk.x,      control = ctrl)
covMcd(log(brain), control = ctrl)
}
\keyword{robust}
\keyword{multivariate}
