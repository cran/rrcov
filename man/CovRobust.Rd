\name{CovRobust}
\alias{CovRobust}
\title{ Robust Location and Scatter Estimation }
\description{
    Computes a robust multivariate location and scatter estimate with a high
     breakdown point, using one of the available estimators.  
}
\usage{
CovRobust(x, control, na.action = na.fail)
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{control}{ a control object (S4) for one of the available control classes, 
    e.g. \code{\link{CovControlMcd-class}}, \code{\link{CovControlOgk-class}}, 
    \code{\link{CovControlSest-class}}, etc., 
    containing estimation options. The class of this object defines 
    which estimator will be used. Alternatively a character string can be specified
    which names the estimator - one of auto, mcd, ogk, m, mve, sfast, surreal, 
    bisquare, rocke. If 'auto' is specified or the argument is missing, the 
    function will select the estimator (see below for details)}
  \item{na.action}{ A function to specify the action to be taken if 
    missing values are found. The default action is for the procedure to fail. 
    An alternative is na.omit, which leads to rejection of cases with missing 
    values on any required variable.}
}
\details{
  This function simply calls the \code{restimate} method of the control object 
  \code{control}. If a character string naming an estimator is specified, a 
  new control object will be created and used. If this argument is missing or 
  'auto' is specified, the function will select the robust estimator 
  (see also \code{\link[robust]{covRob}}). If there are less than 1000 
  observations and less than 10 variables or less than 5000 observations 
  and less than 5 variables, S estimates computed by the \code{sfast} will be used
  (since the Donoho-Stahel estimator is not yet implemented). If there are less than 
  50000 observations and less than 20 variables then the MCD is used. 
  Otherwise the Orthogonalized Quadrant Correlation estimator 
  (CovOgk with the corresponding parameters) is used. 
}
\value{
    An object derived from a \code{CovRobust} object, depending on the selected estimator.
}
%\references{}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovRobust(hbk.x)
CovRobust(hbk.x, CovControlSest(method="bisquare"))

}
\keyword{robust}
\keyword{multivariate}
