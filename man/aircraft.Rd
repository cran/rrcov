\name{aircraft}
\alias{aircraft}
\alias{aircraft.x}
\alias{aircraft.y}
\docType{data}
\title{ Aircraft Data
}
\description{
   Aircraft Data, deals with 23 single-engine aircraft built 
over the years 1947-1979, from Office of Naval Research. 
The dependent variable is cost  (in units of \$100,000) 
and the explanatory variables are aspect ratio, lift-to-drag ratio, 
weight of plane (in pounds) and maximal thrust.

}
\usage{data(aircraft)}
\format{
  A data frame with 23 observations on the following 5 variables.
  \describe{
    \item{\code{X1}}{Aspect Ratio}
    \item{\code{X2}}{Lift-to-Drag Ratio}
    \item{\code{X3}}{Weight}
    \item{\code{X4}}{Thrust}
    \item{\code{Y}}{Cost}
  }
 For convenience, the data sets \code{aircraft.x}, a matrix with 
the four (independent) variables of the data frame, and \code{aircraft.y}, 
the numeric vector giving the fifth (dependent) variable, are provided as well. 
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, page 154, table 22. 

}
\examples{
data(aircraft)
covMcd(aircraft.x)

summary(lm.aircraft <- lm(aircraft.y ~ aircraft.x))

}
\keyword{datasets}
