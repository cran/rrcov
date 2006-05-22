\name{CovRobust-class}
\docType{class}
\alias{CovRobust-class}
\alias{isClassic,CovRobust-method}
\alias{show,CovRobust-method}
\alias{summary,CovRobust-method}

\title{Class "CovRobust" - virtual base class for robust
 estimates of multivariate location and scatter}
\description{  \code{CovRobust} is a virtual base class used for deriving the concrete classes
representing different robust estimates of multivariate location and scatter. Here are implemeted the
standard methods common for all robust estimates like \code{show}, \code{summary} and \code{plot}.
The derived classes can override these methods and can define new ones.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
    \item{\code{center}:}{location}
    \item{\code{cov}:}{covariance matrix}
    \item{\code{n.obs}:}{number of observations used to compute the estimates}
    \item{\code{iter}:}{number of iterations used to compute the estimates}
    \item{\code{crit}:}{value of the criterion function}
    \item{\code{wt}:}{weights}
    \item{\code{call}:}{the call to the function that returns the object}
    \item{\code{corr}:}{whether the correlation matrix was computed}
    \item{\code{mah}:}{robust distances}
    \item{\code{method}:}{a character string describing the method used to compute the estimate.}
    \item{\code{X}:}{data}
  }
}
\section{Extends}{
    Class \code{"Cov"}, directly.
}
\section{Methods}{
\describe{
\item{isClassic}{\code{signature(obj = "CovRobust")}: Will return FALSE, since this is a 'Robust' object }
\item{show}{\code{signature(object = "CovRobust")}: display the object }
\item{plot}{\code{signature(x = "CovRobust")}: plot the object }
}
}
%\references{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}

\seealso{
  \code{\link{Cov}}, \code{\link{Cov-class}}, \code{\link{CovMest-class}}
}
\examples{
     data(hbk)
     hbk.x <- data.matrix(hbk[, 1:3])
     cv <- CovMest(hbk.x)               # it is not possible to create an object of 
                                        # class CovRobust, since it is a VIRTUAL class
     cv
     summary(cv)                        # summary method for class CovRobust
     plot(cv)                           # plot method for class CovRobust
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
