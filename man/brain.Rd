\name{brain}
\alias{brain}
\docType{data}
\title{ Brain and Body Weights for 28 Species }
\description{
   Average brain and body weights for 28 species of land animals. 
If the log of the two variables, i.e. log(brain weight) and 
log(body weight) are taken, the classical Mahalanobis distance 
declares only observation 25 (Brachiosaurus) as an outlier, 
whereas the robust distances computed using for example MCD - 
covMcd() - declare four observations as outliers (25, 6, 14 and 16). 
These are three dinosaurs with relatively small brains and the 
human with relatively heavy brain.
}
\usage{data(brain)}
\format{
  A data frame with 28 observations on the following 2 variables.
  \describe{
    \item{\code{body}}{body weight in kg}
    \item{\code{brain}}{brain weight in g}
  }
}
\Note{
  This data set is also available in the package MASS under the name 'Animals', but the observations 25 (Brachiosaurus) and 26 (Rat) are swaped.
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.57. 
 
}
\examples{
data(brain)
lbrain<-log(brain)
plot(mahalanobis(lbrain,mean(lbrain),var(lbrain)))

mcd<-covMcd(lbrain)
plot(mahalanobis(lbrain,mcd$center,mcd$cov))

}
\keyword{datasets}
