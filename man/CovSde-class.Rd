\name{CovSde-class}
\Rdversion{1.1}
\docType{class}
\alias{CovSde-class}

\title{
Stahel-Donoho Estimates of Multivariate Location and Scatter
}
\description{ This class, derived from the virtual class \code{"CovRobust"} 
    accomodates Stahel-Donoho estimates of multivariate location and scatter.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovSde", ...)},
but the usual way of creating \code{CovSde} objects is a call to the function
\code{CovSde} which serves as a constructor.
}
\section{Slots}{
   \describe{
    \item{\code{center}:}{Object of class \code{"vector"} - 
        the estimate of the location vector}
    \item{\code{cov}:}{Object of class \code{"matrix"} - 
        the estimate of the covariance matrix}
    \item{\code{n.obs}:}{Object of class \code{"numeric"} - number of observations 
        used for the estimates}
    \item{\code{mah}:}{Object of class \code{"Uvector"} - vector of distances 
        computed relative to \code{location} and \code{center} or NULL}
    \item{\code{wt}:}{Object of class \code{"vector"} - vector of weights (w)}
    \item{\code{iter}:}{Object of class \code{"numeric"} - number of iterations }
    \item{\code{crit}:}{Object of class \code{"numeric"} - the value of the objective function}
    \item{\code{call}:}{Object of class \code{"language"} - the call to \code{CovMest}}
    \item{\code{method}:}{Object of class \code{"character"} - method used = \code{"Stahel-Donoh Estimates"}}
    \item{\code{singularity}:}{a list with singularity information for the 
        covariance matrix (or \code{NULL} if not singular)}
    \item{\code{X}:}{Object of class \code{"Umatrix"} - the input data matrix or NULL }
    
  }
}
\section{Extends}{
Class \code{"\linkS4class{CovRobust}"}, directly.
Class \code{"\linkS4class{Cov}"}, by class "CovRobust", distance 2.
}
\section{Methods}{
No methods defined with class "CovSde" in the signature.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovSde}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovSde")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
