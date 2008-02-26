\name{CovMve}
\alias{CovMve}
%
\title{ Robust Location and Scatter Estimation via MVE }
\description{
     Computes a robust multivariate location and scatter estimate with a high
     breakdown point, using the \sQuote{Fast MVE} (Minimum Volume Ellipsoid) 
     estimator.
}
\usage{
CovMve(x, alpha = 1/2, nsamp = 500, seed = NULL, trace = FALSE, control)
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{nsamp}{ number of subsets used for initial estimates or \code{"best"}
    or \code{"exact"}.  Default is \code{nsamp = 500}.  For
    \code{nsamp="best"} exhaustive enumeration is done, as long as the
    number of trials does not exceed 5000. For \code{"exact"},
    exhaustive enumeration will be attempted however many samples are
    needed.  In this case a warning message will be displayed saying
    that the computation can take a very long time.}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{control}{ a control object (S4) of class \code{\link{CovControlMve-class}} 
    containing estimation options - same as these provided in the fucntion 
    specification. If the control object is supplied, the parameters from it 
    will be used. If parameters are passed also in the invocation statement, they will 
    override the corresponding elements of the control object.}
}
\details{
  This function computes the minimum covariance determinant estimator 
  of location and scatter and returns an S4 object of class 
  \code{\link{CovMve-class}} containing the estimates. 
}
\value{
  An S4 object of class \code{\link{CovMve-class}} which is a subclass of the 
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}
\seealso{
  \code{\link[MASS]{cov.mve}} from package \pkg{MASS}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovMve(hbk.x)

## the following three statements are equivalent
c1 <- CovMve(hbk.x, alpha = 0.75)
c2 <- CovMve(hbk.x, control = CovControlMve(alpha = 0.75))
## direct specification overrides control one:
c3 <- CovMve(hbk.x, alpha = 0.75,
             control = CovControlMve(alpha=0.95))
c1
}
\keyword{robust}
\keyword{multivariate}
