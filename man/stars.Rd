\name{stars}
\alias{stars}
\alias{stars.x}
\alias{stars.y}
\docType{data}
\title{ Data for the Hertzsprung-Russell Diagram of the star cluster CYG OB1 }
\description{
   Data for the Hertzsprung-Russell Diagram of the Star Cluster CYG OB1, 
which contains 47 stars in the direction of Czgnus, from C.Doom.
The first variable is the logarithm of the effective temperature at the surface of 
the star (Te) and the second one is the logarithm of its light intencity (L/L0).

In the Hertzsprung-Russell diagram, which is the scatterplot of these data points , 
where the log temperature is ploted from left to right, two groups of points are seen: 
the majority which tend to follow a steep band and four stars in the upper corner. 
In the astronomy the 43 stars are said to lie on the main sequence and the four 
remaining stars are called giants (the points 11, 20, 30, 34).

}
\usage{data(stars)}
\format{
  A data frame with 47 observations on the following 2 variables.
  \describe{
    \item{\code{X1}}{Logarithm of the effective temperature at the surface of the star (Te)}
    \item{\code{Y}}{Logarithm of its light intencity (L/L0)}
  }

 For convenience, the data sets \code{stars.x}, a matrix 
 with the first (independent) variable of the data frame, 
 and \code{stars.y}, the numeric vector giving the second (dependent) 
 variable, are provided as well. 

}
%\Note{}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.27, table 3. 
 
}
\examples{
data(stars)
plot(stars)
covMcd(stars)
lm.stars <- lm(stars.y ~ stars.x)
summary(lm.stars)
plot(lm.stars)
lts.stars<-ltsReg(stars.x, stars.y)
plot(lts.stars)

}
\keyword{datasets}
