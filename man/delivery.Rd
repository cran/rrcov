\name{delivery}
\alias{delivery}
\alias{delivery.x}
\alias{delivery.y}
\docType{data}
\title{ Delivery Time Data, from Montgomery and Peck (1982)
}
\description{
   Delivery Time Data, from Montgomery and Peck (1982)
The aim is to explain the time required to service a vending 
machine by means of the number of products stocked and the 
distance walked by the route driver.
}
\usage{data(delivery)}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{X1}}{Number of Products}
    \item{\code{X2}}{Distance}
    \item{\code{Y}}{Delivery time }
  }
 For convenience, the data sets \code{delivery.x}, a matrix with 
the two (independent) variables of the data frame, and \code{delivery.y}, 
the numeric vector giving the thirth (dependent) variable, are provided as well. 
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, page 155, table 23. 

}
\examples{
data(delivery)
covMcd(delivery.x)

summary(lm.delivery <- lm(delivery.y ~ delivery.x))

}
\keyword{datasets}
