\name{phosphor}
\alias{phosphor}
\alias{phosphor.y}
\alias{phosphor.x} 
\docType{data}
\title{ Phosphorus Content Data }
\description{
   This dataset investigates the effect from inorganic and organic Phosphorus in the soil upon the phosphorus content of the corn grown in this soil, from Prescott (1975).
}
\usage{data(phosphor)}
\format{
  A data frame with 18 observations on the following 3 variables.
  \describe{
    \item{X1}{Inorganic Phosphorus }
    \item{X2}{Organic Phosphorus }
    \item{Y}{Plant Phosphorus }
  }
  
 For convenience, the data sets \code{phosphor.x}, a matrix with the first two (independent) variables of the data frame, and \code{phosphor.y}, the numeric vector giving the thirth (dependent) variable, are provided as well. 

}
%\Note{} 
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.156, table 24. 
 
}

%\details{} 

\examples{
data(phosphor)
plot(phosphor.x)
covMcd(phosphor.x)
summary(lm.phosphor <- lm(phosphor.y ~ phosphor.x))
}
\keyword{datasets}
