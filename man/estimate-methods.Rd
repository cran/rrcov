\name{estimate-methods}
\docType{methods}
\alias{estimate}
\alias{estimate-methods}
\alias{estimate,CovControlMest-method}
\title{ Methods for Function estimate in Package 'rr4cov'}
\description{
    Each concrete control class, like \code{CovControlMest}, \code{CovControlOkg}, 
    etc., should implement an \code{estimate} method which will call the correponding
    (constructor)-function and will return the obtained S4 class, derived from
    \code{CovRobust}.
}
\usage{
\S4method{estimate}{CovControlMest}(obj, x) 
}
\arguments{
 \item{obj}{an object of class \code{"CovControlEstimate"}}
  \item{x}{ Data frame or matrix containing the data }.
  \item{... }{ other parameters to be passed through to the estimation function. }
}


\section{Methods}{
\describe{

\item{obj = "CovControlMest"}{ Compute the constrained M-estimates by calling\code{CovMest} }
}}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
