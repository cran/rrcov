\name{CovSest-class}
\docType{class}
\alias{CovSest-class}

\title{ S Estimates of Multivariate Location and Scatter }
\description{ This class, derived from the virtual class \code{"CovRobust"} 
    accomodates S Estimates of multivariate location and scatter computed 
    by the  \sQuote{Fast S} or \sQuote{SURREAL} algorithm.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovSest", ...)},
but the usual way of creating \code{CovSest} objects is a call to the function
\code{CovSest} which serves as a constructor.
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
    \item{\code{method}:}{Object of class \code{"character"} - method used = \code{"S-Estimates"}}
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
No methods defined with class "CovSest" in the signature.
}
%\references{}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovSest}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovSest")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
