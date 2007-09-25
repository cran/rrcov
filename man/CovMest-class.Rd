\name{CovMest-class}
\docType{class}
\alias{CovMest-class}

\title{Constrained M-estimates of Multivariate Location and Scatter }
\description{ 
    This class, derived from the virtual class "CovRobust" accomodates
    constrained M-Estimates of multivariate location and scatter 
    based on the translated biweight function (\sQuote{t-biweight}) using
    a High breakdown point initial estimate (Minimum Covariance Determinant - 
    \sQuote{Fast MCD})}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovMest", ...)},
    but the usual way of creating \code{CovMest} objects is a call to the function
    \code{CovMest} which serves as a constructor.

}
\section{Slots}{
     \describe{
    \item{\code{center}:}{Object of class \code{"vector"} - the estimate of the location vector}
    \item{\code{cov}:}{Object of class \code{"matrix"} - the estimate of the covariance matrix}
    \item{\code{n.obs}:}{Object of class \code{"numeric"} - number of observations used for the estimates}
    \item{\code{mah}:}{Object of class \code{"Uvector"} - vector of distances computed relative to \code{location} and \code{center} or NULL}
    \item{\code{iter}:}{Object of class \code{"numeric"} - number of iterations }
    \item{\code{wt}:}{Object of class \code{"vector"} - vector of weights (w)}
    \item{\code{vt}:}{Object of class \code{"vector"} - vector of weights (v)}
    \item{\code{crit}:}{Object of class \code{"numeric"} - the value of the objective function}
    \item{\code{call}:}{Object of class \code{"language"} - the call to \code{CovMest}}
    \item{\code{method}:}{Object of class \code{"character"} - method used = \code{"M-Estimates"}}
    \item{\code{singularity}:}{a list with singularity information for the 
        ocvariance matrix (or \code{NULL} of not singular)}
    \item{\code{X}:}{Object of class \code{"Umatrix"} - the data matrix or NULL }
  }
}
\section{Extends}{
   Class \code{"CovRobust"}, directly.
   Class \code{"Cov"}, by class \code{"CovRobust"}.
}
\section{Methods}{
   No methods defined with class "CovMest" in the signature.
}
\references{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\note{}

\seealso{
   \code{\link{CovMest}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
