\name{brain}
\alias{brain}
\docType{data}
\title{ Brain and Body Weights for 28 Species }
\description{
   Average brain and body weights for 28 species of land animals.
}
\usage{data(brain)}
\format{
  A data frame with 28 observations on the following 2 variables.
  \describe{
    \item{body}{body weight in kg}
    \item{brain}{brain weight in g}
  }
}
\Note{
  This data set is also available in the package MASS as 'Animals'
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley, p.57. 
 
}
\examples{
data(brain)
plot(brain)
covMcd(brain)
}
\keyword{datasets}
