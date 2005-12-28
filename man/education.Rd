\name{education}
\alias{education}
\alias{education.x}
\alias{education.y}
\docType{data}
\title{ Education Expenditure Data
 }
\description{
   Education Expenditure Data, from Chatterjee and Price (1977, p.108). This data set, representing the education expenditure variables in the 50 US states,  provides an interesting example of heteroscedacity. 
}
\usage{data(education)}
\format{
  A data frame with 50 observations on the following 6 variables.

  \describe{
    \item{\code{State}}{State}
    \item{\code{Region}}{Region (1=Northeastern, 2=North central, 3=Southern, 4=Western)}
    \item{\code{X1}}{Number of residents per thousand residing in urban areas in 1970}
    \item{\code{X2}}{Per capita personal income in 1973}
    \item{\code{X3}}{Number of residents per thousand under 18 years of age in 1974}
    \item{\code{Y}}{Per capita expenditure on public education in a state, projected for 1975}
  }

 For convenience, the data sets \code{education.x}, a matrix with the 
 three (explanatory) variables of the data frame, and \code{education.y}, 
 the numeric vector giving the dependent variable, are provided as well. 

}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.110, table 16. 

}
\examples{
data(education)
covMcd(education.x)
summary(lm.education <- lm(Y~X1:X3, data=education))
summary(lts.education <- ltsReg(Y~X1:X3, data=education))
}
\keyword{datasets}
