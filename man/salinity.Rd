\name{salinity}
\alias{salinity}
\alias{salinity.x}
\alias{salinity.y}
\docType{data}
\title{ Salinity Data
 }
\description{
   That is a set of measurements of water salinity (i.e., its salt concentration) and river discharge in taken in North Carolina's Pamlico Sound; from Ruppert and Carroll (1980).

}
\usage{data(salinity)}
\format{
  A data frame with 28 observations on the following 4 variables.
  \describe{
    \item{X1}{Lagged Salinity}
    \item{X2}{Trend}
    \item{X3}{Discharge}
    \item{Y}{Salinity}
  }

For convenience, the data sets \code{salinity.x}, a matrix with the three (independent) variables of the data frame, and \code{salinity.y}, the numeric vector giving the fourth (dependent) variable, are provided as well. 
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.57. 

}
\examples{
data(salinity)
covMcd(salinity.x)

summary(lm.salinity <- lm(salinity.y ~ salinity.x))

}
\keyword{datasets}
