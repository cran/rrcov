\name{CovMcd}
\alias{CovMcd}
%
\title{ Robust Location and Scatter Estimation via MCD }
\description{
     Computes a robust multivariate location and scatter estimate with a high
     breakdown point, using the \sQuote{Fast MCD} (Minimum Covariance
     Determinant) estimator.
}
\usage{
CovMcd(x, alpha = 1/2, nsamp = 500, seed = NULL, trace = FALSE, use.correction = TRUE, control)
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
  \item{use.correction}{ whether to use finite sample correction factors.
    Default is \code{use.correction=TRUE}}
  \item{control}{ a control object (S4) of class \code{\link{CovControlMcd-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
}
\details{
  This function computes the minimum covariance determinant estimator
  of location and scatter and returns an S4 object of class
  \code{\link{CovMcd-class}} containing the estimates.
  The implementation of the function is similar to the existing \R function
  \code{\link[robustbase]{covMcd}()} which returns an S3 object.
  The MCD method looks for the \eqn{h (> n/2)}
  observations (out of \eqn{n}) whose classical
  covariance matrix has the lowest possible determinant.  The raw MCD
  estimate of location is then the average of these \eqn{h} points,
  whereas the raw MCD estimate of scatter is their covariance matrix,
  multiplied by a consistency factor and a finite sample correction factor
  (to make it consistent at the normal model and unbiased at small samples).
  Both rescaling factors are returned also in the vector \code{raw.cnp2}
  of length 2. Based on these raw MCD estimates, a reweighting step is performed
  which increases the finite-sample efficiency considerably - see Pison et al. (2002).
  The rescaling factors for the reweighted estimates are returned in the
  vector \code{cnp2} of length 2. Details for the computation of the finite
  sample correction factors can be found in Pison et al. (2002).
  The finite sample corrections can be suppressed by setting \code{use.correction=FALSE}.
  The implementation in rrcov uses the Fast MCD algorithm of Rousseeuw and Van Driessen (1999)
  to approximate the minimum covariance determinant estimator.
}
\value{
  An S4 object of class \code{\link{CovMcd-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.

  Pison, G., Van Aelst, S., and Willems, G. (2002),
  Small Sample Corrections for LTS and MCD,
  \emph{Metrika}, \bold{55}, 111-123.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}
%\note{
%    The function \code{\link[robustbase]{covMcd}()} on which \code{CovMcd} is based
%    will be deprecated i.e. is going to be made obsolete in future versions.
%}
\seealso{
  \code{\link[MASS]{cov.mcd}} from package \pkg{MASS}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovMcd(hbk.x)

## the following three statements are equivalent
c1 <- CovMcd(hbk.x, alpha = 0.75)
c2 <- CovMcd(hbk.x, control = CovControlMcd(alpha = 0.75))
## direct specification overrides control one:
c3 <- CovMcd(hbk.x, alpha = 0.75,
             control = CovControlMcd(alpha=0.95))
c1
}
\keyword{robust}
\keyword{multivariate}
