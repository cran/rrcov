\name{salinity}
\alias{salinity}
\alias{salinity.x}
\alias{salinity.y}
\docType{data}
\title{ Salinity Data
 }
\description{
   This is a data set consisting of measurements of water salinity (i.e., its salt
concentration) and river discharge taken in North Carolina's
Pamlico Sound; This dataset was listed by Ruppert and Carroll
(1980). In Carrol and Ruppert (1985) the physical background of
the data is described. They indicated that observations 5 and 16
correspond to periods of very heavy discharge and showed that the
discrepant observation 5 was masked by observations 3 and 16 i.e.
only after deletion of these observations it was possible to
identify the influential observation 5.

 This data set is a prime example of the masking effect.

}
\usage{data(salinity)}
\format{
  A data frame with 28 observations on the following 4 variables.
  \describe{
    \item{\code{X1}}{Lagged Salinity}
    \item{\code{X2}}{Trend}
    \item{\code{X3}}{Discharge}
    \item{\code{Y}}{Salinity}
  }

For convenience, the data sets \code{salinity.x}, a matrix with the three (independent) variables of the data frame, and \code{salinity.y}, the numeric vector giving the fourth (dependent) variable, are provided as well. 
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.82, table 5. 

}
\examples{
data(salinity)
covMcd(salinity.x)

summary(lm.salinity <- lm(salinity.y ~ salinity.x))

}
\keyword{datasets}
