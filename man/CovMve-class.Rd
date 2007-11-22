\name{CovMve-class}
\docType{class}
\alias{CovMve-class}

\title{ MVE Estimates of Multivariate Location and Scatter }
\description{ 
    This class, derived from the virtual class \code{"CovRobust"} accomodates
    MVE Estimates of multivariate location and scatter computed by the  
    \sQuote{Fast MVE} algorithm.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovMve", ...)},
but the usual way of creating \code{CovMve} objects is a call to the function
\code{CovMve} which serves as a constructor.
}
\section{Slots}{
     \describe{

    \item{\code{center}:}{Object of class \code{"vector"} - the final estimate of the location vector}
    \item{\code{cov}:}{Object of class \code{"matrix"} - the final estimate of the covariance matrix}
    \item{\code{n.obs}:}{Object of class \code{"numeric"} - number of observations used for the estimates}
    \item{\code{mah}:}{Object of class \code{"Uvector"} - vector of distances computed 
        relative to the final estimates \code{center} and \code{cov} or NULL}
    \item{\code{iter}:}{Object of class \code{"numeric"} - number of iterations }
    \item{\code{wt}:}{Object of class \code{"vector"} - weights of the observations using the 
        final estimate of the location and scatter}
    \item{\code{crit}:}{Object of class \code{"numeric"} - the value of the objective function - i.e. the determinant}
    \item{\code{call}:}{Object of class \code{"language"} - the call to \code{CovMcd}}
    \item{\code{method}:}{Object of class \code{"character"} - method 
        used = \code{"Minimum Volume Ellipsoid Estimator"}}
    \item{\code{singularity}:}{a list with singularity information for the 
        ocvariance matrix (or \code{NULL} if not singular)}
    \item{\code{alpha}:}{Object of class \code{"numeric"} - the size of the 
        subsets over which the volume of the ellipsoid is minimized (the default is (n+p+1)/2) }
    \item{\code{quan}:}{Object of class \code{"numeric"} - the number of 
        observations on which the MVE is based. If \code{quan} equals 
        \code{n.obs}, the MVE is the classical covariance matrix. }
    \item{\code{best}:}{Object of class \code{"Uvector"} - the best subset 
        found and used for computing the raw estimates. The size of \code{best} 
        is equal to \code{quan} }
    \item{\code{raw.cov}:}{Object of class \code{"matrix"} the raw 
        (not reweighted) estimate of location }
    \item{\code{raw.center}:}{Object of class \code{"vector"} - the raw 
        (not reweighted) estimate of scatter}
    \item{\code{raw.mah}:}{Object of class \code{"Uvector"} - mahalanobis 
        distances of the observations based on the raw estimate of the 
        location and scatter }
    \item{\code{raw.wt}:}{Object of class \code{"Uvector"} - weights of 
        the observations based on the raw estimate of the location and scatter }
    \item{\code{X}:}{Object of class \code{"Umatrix"} - the data matrix or NULL }  }
}
\section{Extends}{
Class \code{"CovRobust"}, directly.
Class \code{"Cov"}, by class \code{"CovRobust"}.
}
\section{Methods}{
No methods defined with class "CovMve" in the signature.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovMve}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovMve")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
