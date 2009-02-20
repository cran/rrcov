\name{CovMcd-class}
\docType{class}
\alias{CovMcd-class}

\title{ MCD Estimates of Multivariate Location and Scatter }
\description{ 
    This class, derived from the virtual class \code{"CovRobust"} accomodates
    MCD Estimates of multivariate location and scatter computed by the  
    \sQuote{Fast MCD} algorithm.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovMcd", ...)},
but the usual way of creating \code{CovMcd} objects is a call to the function
\code{CovMcd} which serves as a constructor.
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
        used = \code{"Minimum Covariance Determinant Estimator"}}
    \item{\code{singularity}:}{a list with singularity information for the 
        covariance matrix (or \code{NULL} if not singular)}
    \item{\code{alpha}:}{Object of class \code{"numeric"} - the size of the 
        subsets over which the determinant is minimized (the default is (n+p+1)/2) }
    \item{\code{quan}:}{Object of class \code{"numeric"} - the number of 
        observations on which the MCD is based. If \code{quan} equals 
        \code{n.obs}, the MCD is the classical covariance matrix. }
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
    \item{\code{raw.cnp2}:}{Object of class \code{"numeric"} - a vector of length 
        two containing the consistency correction factor and the finite sample 
        correction factor of the raw estimate of the covariance matrix }
    \item{\code{cnp2}:}{Object of class \code{"numeric"} - a vector of length two 
        containing the consistency correction factor and the finite sample 
        correction factor of the final estimate of the covariance matrix.
    }
    \item{\code{X}:}{Object of class \code{"Umatrix"} - the data matrix or NULL }
  }
}
\section{Extends}{
Class \code{"CovRobust"}, directly.
Class \code{"Cov"}, by class \code{"CovRobust"}.
}
\section{Methods}{
No methods defined with class "CovMcd" in the signature.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovMcd}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovMcd")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
