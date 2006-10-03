\name{CovOgk-class}
\docType{class}
\alias{CovOgk-class}

\title{ OGK Estimates of Multivariate Location and Scatter }
\description{
    This class, derived from the virtual class \code{"CovRobust"} accomodates
    OGK Estimates of multivariate location and scatter computed by the  
    algorithm proposed by Marona and Zamar (2002).
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovOgk", ...)},
    but the usual way of creating \code{CovOgk} objects is a call to the function
    \code{CovOgk} which serves as a constructor.
}
\section{Slots}{
    \describe{
    \item{\code{center}:}{Object of class \code{"vector"} - the final (reweighted) estimate of the location vector}
    \item{\code{cov}:}{Object of class \code{"matrix"} - the final (reweighted) estimate of the covariance matrix}
    \item{\code{n.obs}:}{Object of class \code{"numeric"} - number of observations used for the estimates}
    \item{\code{mah}:}{Object of class \code{"Uvector"} - vector of distances computed 
        relative to the final estimates \code{center} and \code{cov} or NULL. Use the accessor method 
        \code{getDistance()}.}
    \item{\code{wt}:}{Object of class \code{"vector"} - weights of the observations using the 
        final estimate of the location and scatter (currently not used) }        
    \item{\code{iter}:}{Object of class \code{"numeric"} - number of ortogonalization iterations used}
    \item{\code{crit}:}{not used}
    \item{\code{call}:}{Object of class \code{"language"} - the call to \code{CovOgk}}
    \item{\code{method}:}{Object of class \code{"character"} - method 
        used = \code{"Orthogonalized Gnanadesikan-Kettenring Estimator"}}

    \item{\code{raw.cov}:}{Object of class \code{"matrix"} the raw 
        (not reweighted) estimate of covariance matrix }
    \item{\code{raw.center}:}{Object of class \code{"vector"} - the raw 
        (not reweighted) estimate of the location vector}
    \item{\code{raw.mah}:}{Object of class \code{"Uvector"} - mahalanobis 
        distances of the observations based on the raw estimate of the 
        location and scatter }
    \item{\code{raw.wt}:}{Object of class \code{"Uvector"} - weights of 
        the observations based on the raw estimate of the location and scatter }
    \item{\code{X}:}{Object of class \code{"Umatrix"} - the data matrix or NULL }
  }
}
\section{Extends}{
Class \code{"CovRobust"}, directly.
Class \code{"Cov"}, by class \code{"CovRobust"}.
}
\section{Methods}{
No methods defined with class "CovOgk" in the signature.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}
%\note{}
\seealso{ \code{\link{CovMcd-class}}, \code{\link{CovMest-class}} }
%\examples{}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
