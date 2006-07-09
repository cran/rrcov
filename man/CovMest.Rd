\name{CovMest}
\alias{CovMest}
\title{Constrained M-Estimates of Location and Scatter}
\concept{High breakdown point}
\usage{
    CovMest(x, r = 0.45, arp = 0.05, eps=1e-3, 
        maxiter=120, control, t0, S0)
}
\description{
    Computes constrained M-Estimates of multivariate location and scatter 
    based on the translated biweight function (\sQuote{t-biweight}) using
    a High breakdown point initial estimate (Minimum Covariance Determinant - 
    \sQuote{Fast MCD})
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{r}{ required breakdown point.  Allowed values are between 
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45}
  \item{arp}{ asympthotic rejection point, i.e. the fraction of points 
    receiving zero weight (see Rocke (1996)).  Default is \code{0.05}.}
  \item{eps}{ a numeric value specifying the relative precision of the solution of 
    the M-estimate. Defaults to \code{1e-3}}
  \item{maxiter}{ maximum number of iterations allowed in the computation of 
    the M-estimate. Defaults to 120}
  \item{control}{ a control object (S4) of class \code{\link{CovControlMest-class}} 
    containing estimation options - same as these provided in the fucntion 
    specification. If the control object is supplied, the parameters from it 
    will be used. If parameters are passed also in the invocation statement, they will 
    override the corresponding elements of the control object.}
  \item{t0}{ optional initial high breakdown point estimates of the location. 
    If not supplied MCD will be used. }
  \item{S0}{ optional initial high breakdown point estimates of the 
    scatter. If not supplied MCD will be used. }
}
\details{
  Rocke (1996) has shown that the S-estimates of multivariate location and scatter
  in high dimensions can be sensitive to outliers even if the breakdown point 
  is set to be near 0.5. To mitigate this problem he proposed to utilize 
  the translated  biweight (or t-biweight) method with a
  standardization step consisting of equating the median of \code{rho(d)} 
  with the median under normality. This is then not an S-estimate, but is 
  instead a constrained M-estimate. In order to make the smooth estimators 
  to work, a reasonable starting point is necessary, which will lead reliably to a 
  good solution of the estimator. In \code{covMest} the MCD computed by 
  \code{\link{covMcd}} is used, but the user has the possibility to give her own 
  initial estimates.
  
}
\note{
    The psi, rho and weight functions for the M estimation are encapsulated in a 
    virtual S4 class \code{PsiFun} from which a \code{PsiBwt} class, implementing 
    the translated biweight (t-biweight), is dervied. The base class  \code{PsiFun} 
    contains also the M-iteration itself. Although not documented and not 
    accessibale directly by the user these classes will form the bases for adding 
    other functions (biweight, LWS, etc.) as well as S-estimates.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}, 

    (some code from C. Becker -  
    http://www.sfb475.uni-dortmund.de/dienst/de/content/struk-d/bereicha-d/tpa1softw-d.html)
}
\value{
  An object of class \code{CovMest} which is a subclass of the virtual class \code{CovRobust}.
}
\references{
    D.L.Woodruff and D.M.Rocke (1994) Computable robust estimation of multivariate location
    and shape on high dimension using compound estimators, \emph{Journal of the American 
    Statistical Association}, \bold{89}, 888--896.
    
    D.M.Rocke (1996) Robustness properties of S-estimates of multivariate location and
    shape in high dimension, \emph{Annals of Statistics}, \bold{24}, 1327-1345.
    
    D.M.Rocke and D.L.Woodruff (1996) Identification of outliers in multivariate data \emph{Journal
    of the American Statistical Association}, \bold{91}, 1047--1061.
}
\seealso{
  \code{\link{covMcd}}, \code{\link{Cov}}, \code{\link{Cov-class}}, 
  \code{\link{CovRobust-class}}, \code{\link{CovMest-class}}
}
\examples{

library(rrcov)
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovMest(hbk.x)

## the following four statements are equivalent
c0 <- CovMest(hbk.x)
c1 <- CovMest(hbk.x, r = 0.45)
c2 <- CovMest(hbk.x, control = CovControlMest(r = 0.45))
c3 <- CovMest(hbk.x, control = new("CovControlMest", r = 0.45))

## direct specification overrides control one:
c4 <- CovMest(hbk.x, r = 0.40,
             control = CovControlMest(r = 0.25))
c1
summary(c1)
plot(c1)
}
\keyword{robust}
\keyword{multivariate}