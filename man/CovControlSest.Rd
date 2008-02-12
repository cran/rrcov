\name{CovControlSest}
\alias{CovControlSest}
\title{ Constructor function for objects of class "CovControlSest" }
\concept{High breakdown point}
\description{
    This function will create a control object \code{CovControlSest} 
    containing the control parameters for \code{CovSest}
}
\usage{
    CovControlSest(bdp = 0.5, nsamp = 500, seed = NULL, trace = FALSE, tolSolve = 1e-13, algo = "sfast")
}
\arguments{
  \item{bdp}{ a numeric value specifying the required 
    breakdown point.  Allowed values are between 
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45 
  }
  \item{nsamp}{the number of random subsets considered. Default is \code{nsamp = 500}.}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}  
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.
  }
  \item{algo}{ Which algorithm to use: 'sfast'=FAST-S or 'surreal'=SURREAL }
}
%\details{}
\value{
A \code{CovControlSest} object
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\seealso{}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlSest", bdp=0.4)
    ctrl2 <- CovControlSest(bdp=0.4)

    data(hbk)
    CovSest(hbk, control=ctrl1)    
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
