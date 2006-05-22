\name{getCenter-methods}
\docType{methods}
\title{Accessor methods to the essential slots of Cov and its subclasses}
\alias{getCenter}
\alias{getCenter-methods} 
\alias{getCov}
\alias{getCov-methods} 
\alias{getCorr}
\alias{getCorr-methods} 
\alias{getData}
\alias{getData-methods} 
\alias{getDistance}
\alias{getDistance-methods} 
\alias{getEvals}
\alias{getEvals-methods} 
\alias{isClassic}
\alias{isClassic-methods} 
\usage{
getCenter(obj)
getCov(obj)
getCorr(obj)
getData(obj)
getDistance(obj)
getEvals(obj)
isClassic(obj)
} 

\arguments{
    \item{obj}{an object of class \code{"Cov"}  or of a class derived from \code{"Cov"}}
}

\description{
    Accessor methods to the slots of objects of class\code{Cov} and its subclasses}
}
\section{Methods}{
\describe{

\item{obj = "Cov"}{generic functions - see \code{getCenter}, \code{getCov}, \code{getCorr}, \code{getData}, \code{getDistance}, \code{getEvals}, \code{isClassic}}
\item{obj = "CovRobust"}{generic functions - see \code{getCenter}, \code{getCov}, \code{getCorr}, \code{getData}, \code{getDistance}, \code{getEvals}, \code{isClassic}}
}}
\keyword{methods}
\keyword{multivariate}
\keyword{robust}
