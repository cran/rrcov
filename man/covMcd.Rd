\name{covMcd}
\alias{covMcd}
\alias{print.mcd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust location and scatter estimation with high breakdown point }
\description{
     Compute a multivariate location and scale estimate with a high
     breakdown point using the Fast MCD (Minimum Covariance Determinant) Estimator.
}
\usage{
covMcd(x, cor=FALSE, alpha=1/2, nsamp=500, seed=0, print.it=FALSE, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a matrix or data frame. }
  \item{cor}{ should the returned result include a correlation matrix? Default is \code{cor = FALSE} }.
  \item{alpha}{ This parameter controls the size of the subsets over which the determinant is
minimized, i.e. \code{alpha*n} observations are used for computing the determinant. Allowed values are between 0.5 and 1 and the default is 0.5. }
  \item{nsamp}{ number of subsets used for initial estimates. Default is \code{nsamp = 500}}
  \item{seed}{ starting value for random generator. Default is \code{seed = 0}}
  \item{print.it}{ whether to print intermediate results. Default is \code{print.it = FALSE}}
  \item{control}{ a list with estimation options - same as these provided in the 
  fucntion specification. If the control object is supplied, the parameters from it 
  will be used. If parameters are passed also in the invocation statement, they will 
  override the corresponding elements of the control object. }
}
\details{
  The minimum covariance determinant estimator of location and
scatter implemented in covMcd() is similar to the existing 
R function cov.mcd() in MASS. The MCD method looks for the
h(> n/2) observations (out of n) whose classical covariance 
matrix has the lowest possible determinant. The raw MCD estimate 
of location is then the average of these h points, whereas the raw MCD 
estimate of scatter is their covariance matrix, multiplied with a
consistency factor. Based on these raw MCD estimates, a reweighting step is performed
which increases the finite-sample eficiency considerably - see Pison et.al. (2002). 
The implementation in rrcov uses the Fast MCD algorithm of Rousseeuw and Van Driessen (1999) to
approximate the minimum covariance determinant estimator. 

}
\value{
  A list with components

  \item{center}{
    the final estimate of location.
  }
  \item{cov}{
    the final estimate of scatter.
  }
  \item{cor}{
    the (final) estimate of the correlation matrix (only if \code{cor = TRUE}) .
  }
  \item{crit}{
    the value of the criterion, i.e. the determinant.
  }
  \item{best}{
    the best subset found and used for computing the raw estimates. The size of \code{best} is equal to \code{quan}.
  }
  \item{mah}{
    mahalanobis distances of the observations using the final estimate of the location and scatter.
  }
  \item{mcd.wt}{
    weights of the observations using the final estimate of the location and scatter.
  }
  \item{raw.center}{
    the raw (not reweighted) estimate of location.
  }
  \item{raw.cov}{
    the raw (not reweighted) estimate of scatter.
  }
  \item{raw.mah}{
    mahalanobis distances of the observations based on the raw estimate of the location and scatter.
  }
  \item{raw.weights}{
    weights of the observations based on the raw estimate of the location and scatter.
  }
  \item{X}{
    the input data as a matrix.
  }
  \item{n.obs}{
    total number of observations. 
  }
  \item{alpha}{ 
    the size of the subsets over which the determinant is minimized (the default is (n+p+1)/2). 
  }  
  \item{quan}{ 
    the number of observations on which the MCD is based. 
    If \code{quan} equals \code{n.obs}, the MCD is the classical covariance matrix.
  }  
  \item{method}{
    character string naming the method (Minimum Covariance Determinant).
  }
  
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

% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{}
\examples{

data(hbk)
covMcd(hbk.x)

# the following three statements are equivalent
covMcd(hbk.x, alpha=0.75)
covMcd(hbk.x, control = rrcov.control(alpha=0.75))
covMcd(hbk.x, alpha = 0.75, control = rrcov.control(alpha=0.95))

}
\keyword{robust}
\keyword{multivariate}
