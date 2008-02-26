\encoding{latin1}
\name{CovSest}
\alias{CovSest}
\title{ S Estimates of Multivariate Location and Scatter }
\concept{High breakdown point}
\description{
    Computes S-Estimates of multivariate location and scatter based on Tukey's 
    biweight function using a fast algorithm similar to the one proposed by 
    Salibian-Barrera and Yohai (2006) for the case of regression. Alternativley, the 
    Ruppert's SURREAL algorithm can be used.
}
\usage{
    CovSest(x, nsamp = 500, bdp = 0.5, seed = NULL, trace = FALSE, tolSolve = 1e-13, 
        algo = c("sfast", "surreal"), control)
}
\arguments{
  \item{x}{ a matrix or data frame. }
  \item{nsamp}{the number of random subsets considered. Default is \code{nsamp = 500}}.
  \item{bdp}{ required breakdown point.  Allowed values are between 
    \code{(n - p)/(2 * n)} and 1 and the default is 0.5}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}  
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.}
  \item{algo}{ Which algorithm to use: 'sfast'=FAST-S or 'surreal'=SURREAL }
  \item{control}{ a control object (S4) of class \code{\link{CovControlSest-class}} 
    containing estimation options - same as these provided in the fucntion 
    specification. If the control object is supplied, the parameters from it 
    will be used. If parameters are passed also in the invocation statement, they will 
    override the corresponding elements of the control object.}
}
\details{
    Computes biweight multivariate S-estimator of location and scatter. The 
    algorithm used is similar to the one proposed by Salibian-Barrera and Yohai (2006) 
    for the case of regression. Alternativley, if the parameter \code{algo} is set to 
    \code{surreal} the Ruppert's SURREAL algorithm will be used.
}
\value{
  An S4 object of class \code{\link{CovSest-class}} which is a subclass of the 
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  H.P. Lopuhaä (1989) On the Relation between S-estimators and M-estimators of 
  Multivariate Location and Covariance.
  \emph{Annals of Statistics} \bold{17} 1662--1683. 

  D. Ruppert (1992) Computing S Estimators for Regression and Multivariate 
  Location/Dispersion. 
  \emph{Journal of Computational and Graphical Statistics} \bold{1} 253--270. 

  M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for {S}-regression
  estimates, \emph{Journal of Computational and Graphical Statistics}, \bold{15},
  414--427.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} and
Matias Salibian-Barrera \email{matias@stat.ubc.ca}.  See also the code from 
Kristel Joossens, K.U. Leuven, Belgium and Ella Roelant, Ghent University, Belgium.
}
%\note{}
%\seealso{}
\examples{

library(rrcov)
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovSest(hbk.x)

## the following four statements are equivalent
c0 <- CovSest(hbk.x)
c1 <- CovSest(hbk.x, bdp = 0.25)
c2 <- CovSest(hbk.x, control = CovControlSest(bdp = 0.25))
c3 <- CovSest(hbk.x, control = new("CovControlSest", bdp = 0.25))

## direct specification overrides control one:
c4 <- CovSest(hbk.x, bdp = 0.40,
             control = CovControlSest(bdp = 0.25))
c1
summary(c1)
plot(c1)

## Use the SURREAL algorithm of Ruppert
cr <- CovSest(hbk.x, algo="surreal")
cr
}
\keyword{robust}
\keyword{multivariate}
