\name{CovOgk}
\alias{CovOgk}
%
\title{ Robust Location and Scatter Estimation - Ortogonalized Gnanadesikan-Kettenring (OGK)}
\description{
     Computes a robust multivariate location and scatter estimate with a high
     breakdown point, using the pairwise algorithm proposed by Marona and 
     Zamar (2002) which in turn is based on the pairwise robust estimator 
     proposed by Gnanadesikan-Kettenring (1972).
}
\usage{
CovOgk(x, niter = 2, beta = 0.9, control)
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{niter}{ number of iterations, usually 1 or 2 since iterations 
    beyond the second do not lead to improvement. }
  \item{beta}{ coverage parameter for the final reweighted estimate }
  \item{control}{ a control object (S4) of class \code{\link{CovControlOgk-class}} 
    containing estimation options - same as these provided in the function 
    specification. If the control object is supplied, the parameters from it 
    will be used. If parameters are passed also in the invocation statement, they will 
    override the corresponding elements of the control object. The control
    object contains also functions for computing the robust univariate location 
    and dispersion estimate \code{mrob} and for computing the robust estimate 
    of the covariance between two random variables \code{vrob}.  }
}
\details{
    The method proposed by Marona and Zamar (2002) allowes to obtain 
    positive-definite and almost affine equivariant robust scatter matrices 
    starting from any pairwise robust scatter matrix. The default robust estimate 
    of covariance between two random vectors used is the one proposed by 
    Gnanadesikan and Kettenring (1972) but the user can choose any other method by 
    redefining the function in slot \code{vrob} of the control object 
    \code{CovControlOgk}. Similarly, the function for computing the robust 
    univariate location and dispersion used is the \code{tau scale} defined 
    in Yohai and Zamar (1998) but it can be redefined in the control object.
    
    The estimates obtained by the OGK method, similarly as in \code{CovMcd} are returned
    as 'raw' estimates. To improve the estimates a reweighting step is performed using
    the coverage parameter \code{beta} and these reweighted estimates are returned as 
    'final' estimates.    
}
\value{
  An S4 object of class \code{\link{CovOgk-class}} which is a subclass of the 
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  Maronna, R.A. and Zamar, R.H. (2002)
  Robust estimates of location and dispersion of high-dimensional datasets;
  \emph{Technometrics} \bold{44}(4), 307--317.

  Yohai, R.A. and Zamar, R.H. (1998)
  High breakdown point estimates of regression by means of the minimization of efficient scale
  \emph{JASA} \bold{86}, 403--413.

  Gnanadesikan, R. and John R. Kettenring (1972)
  Robust estimates, residuals, and outlier detection with multiresponse data.
  \emph{Biometrics} \bold{28}, 81--124.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}
%\note{}
\seealso{ \code{\link{CovMcd}}, \code{\link{CovMest}} }
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovOgk(hbk.x)

## the following three statements are equivalent
c1 <- CovOgk(hbk.x, niter=1)
c2 <- CovOgk(hbk.x, control = CovControlOgk(niter=1))

## direct specification overrides control one:
c3 <- CovOgk(hbk.x, beta=0.95,
             control = CovControlOgk(beta=0.99))
c1

}
\keyword{robust}
\keyword{multivariate}