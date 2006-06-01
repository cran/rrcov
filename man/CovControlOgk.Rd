\name{CovControlOgk}
\alias{CovControlOgk}
%
\title{ Constructor function for objects of class "CovControlOgk" }
\description{
    This function will create a control object \code{CovControlOgk} 
    containing the control parameters for \code{CovOgk}
}
\usage{
CovControlOgk(niter = 2, beta = 0.9, mrob = rrcov:::.mrobTau, vrob = rrcov:::.vrobGK)
}
\arguments{
  \item{niter}{ number of iterations, usually 1 or 2 since iterations 
    beyond the second do not lead to improvement. }
  \item{beta}{ coverage parameter for the final reweighted estimate }
  \item{mrob}{ function for computing the robust univariate location 
    and dispersion - defaults to the \code{tau scale} defined in 
    Yohai and Zamar (1998) 
  }
  \item{vrob}{ function for computing robust estimate 
    of covariance between two random vectors - defaults the one proposed by 
    Gnanadesikan and Kettenring (1972) 
  }
}
%\details{}
\value{
A \code{CovControlMcd} object
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
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\seealso{}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlOgk", beta=0.95)
    ctrl2 <- CovControlOgk(beta=0.95)

    data(hbk)
    CovOgk(hbk, control=ctrl1)    
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
