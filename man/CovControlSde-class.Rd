\name{CovControlSde-class}
\Rdversion{1.1}
\docType{class}
\alias{CovControlSde-class}
\alias{restimate,CovControlSde-method}

\title{Class 'CovControlSde' - contains control parameters for "CovSde"}
\description{ 
    This class extends the \code{CovControl} class 
    and contains the control parameters for CovSde
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovControlSde", ...)} 
    or by calling the constructor-function \code{CovControlSde}.
}
\section{Slots}{
     \describe{
          \item{\code{nsamp}}{a positive integer giving the number of resamples required}
          \item{\code{maxres}}{a positive integer specifying the maximum number of 
            resamples to be performed including those that are discarded due to linearly 
            dependent subsamples.}
          \item{\code{tune = 0.95}}{a numeric value between 0 and 1 giving the fraction of the data to receive non-zero weight.}
          \item{\code{prob = 0.99}}{a numeric value between 0 and 1 specifying the probability of high breakdown point; 
            used to compute \code{nsamp} when \code{nsamp} is omitted.}
          \item{\code{eps = 0.5}}{a numeric value between 0 and 0.5 specifying the breakdown point; 
            used to compute \code{nsamp} when \code{nresamp} is omitted.}
          \item{seed}{starting value for random generator.}
          \item{trace}{whether to print intermediate results.}
          \item{tolSolve}{numeric tolerance to be used for inversion (\code{\link{solve}}) of the covariance matrix in
            \code{\link{mahalanobis}}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{CovControl}"}, directly.
}
\section{Methods}{
  \describe{
    \item{restimate}{\code{signature(obj = "CovControlSde")}: ... }
     }
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlSde", nsamp=2000)
    ctrl2 <- CovControlSde(nsamp=2000)

    data(hbk)
    CovSde(hbk, control=ctrl1)    
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
