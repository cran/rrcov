\name{Cov}
\alias{Cov}
\title{Classical Estimates of Multivariate Location and Scatter}
\concept{High breakdown point}
\description{
    Computes the classical estimates of multivariate location and scatter. 
    Returns an S4 class \code{\link{Cov}} with the estimated \code{center}, 
    \code{cov}, Mahalanobis distances and weights based on these distances. 
    The class \code{Cov} is used as a base class for the classes representing 
    the different robust estimators.
}
\usage{
    Cov(x, unbiased=TRUE)
}
\arguments{
  \item{x}{a matrix or data frame. As usual, rows are observations and
        columns are variables. }
  \item{unbiased}{ whether to return the unbiased estimate of 
    the covariance matrix. Default is  \code{unbiased = TRUE}}
}
\details{
  
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}
}
\value{
  An object of class \code{"Cov"}.   
}
\seealso{
  \code{\link{Cov-class}}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
cv <- Cov(hbk.x)
cv
summary(cv)
plot(cv)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
