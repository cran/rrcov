\name{coleman}
\alias{coleman}
\alias{coleman.x}
\alias{coleman.y}
\docType{data}
\title{ Coleman Data Set
 }
\description{
   Contains information on 20 Schools from the Mid-Atlantic and New England States, 
drawn from a population studied by Coleman et al. (1966). Mosteller and Tukey (1977) 
analyze this sample consisting of measurements on six different variables, one of 
which will be treated as a responce.
}
\usage{data(coleman)}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{\code{X1}}{staff salaries per pupil}
    \item{\code{X2}}{percent of white-collar fathers}
    \item{\code{X3}}{socioeconomic status composite deviation: means for family size, family intactness, father's education, mother's education, and home items}
    \item{\code{X4}}{mean teacher's verbal test score}
    \item{\code{X5}}{mean mother's educational level, one unit is equal to two school years}
    \item{\code{Y}}{verbal mean test score (y, all sixth graders)}
  }

 For convenience, the data sets \code{coleman.x}, a matrix with the 
 five (independent) variables of the data frame, and \code{coleman.y}, 
 the numeric vector giving the sixth (dependent) variable, are provided as well. 

}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.79, table 2. 

}
\examples{
data(coleman)
covMcd(coleman.x)
summary(lm.coleman <- lm(coleman.y ~ coleman.x))

}
\keyword{datasets}
