\name{wood}
\alias{wood}
\alias{wood.x}
\alias{wood.y}
\docType{data}
\title{ Modified Wood Gravity Dataset }
\description{
   
}
\usage{data(wood)}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{X1}{X1}
    \item{X2}{X2}
    \item{X3}{X3}
    \item{X4}{X4}
    \item{X5}{X5}
    \item{Y}{Y}
  }

 For convenience, the data sets \code{wood.x}, a matrix with the first 5 (independent) variables of the data frame, and \code{wood.y}, the numeric vector giving the sixth (dependent) variable, are provided as well. 

}
%\Note{}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley. 
 
}
\examples{
data(wood)
plot(wood.x)
covMcd(wood.x)
summary(lm.wood <- lm(wood.y ~ wood.x))
}
\keyword{datasets}
