\name{heart}
\alias{heart}
\alias{heart.x}
\alias{heart.y}
\docType{data}
\title{ Heart Catherization Data }
\description{
      This data set was analyzed by Weisberg (1980) and Chambers et al. (1983).
   A catheter is passed into a major vein or artery at the femoral region and moved into the heart.
The proper length of the introduced catheter has to be guessed by
the physician. The aim of the data set is to describe the relation
between the catheter length and the patient's height (X1) and
weight (X2).

 This data sets is used to demonstrate the effects caused by collinearity.
The correlation between height and weight is so high that either
variable almost completely determines the other.

}
\usage{data(heart)}
\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{X1}}{Patient's height in inches}
    \item{\code{X2}}{Patient's weights in pounds}
    \item{\code{Y}}{Catheter Length (in centimeters)}
  }

 For convenience, the data sets \code{heart.x}, a matrix 
 with the first two (independent) variables of the data frame, 
 and \code{heart.y}, the numeric vector giving the thirth (dependent) 
 variable, are provided as well. 

}
%\Note{}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.103, table 13. 
 
}
\examples{
data(heart)
plot(heart.x)
covMcd(heart.x)
summary(lm.heart <- lm(heart.y ~ heart.x))
}
\keyword{datasets}
