\name{CovControlMest-class}
\docType{class}
\alias{CovControlMest-class}

\title{Class 'CovControlMest' - contains control parameters for CovMest}
\description{ This class extends the \code{CovControl} class 
    and contains the control parameters for CovMest}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovControlMest", ...)} 
or by calling the constructor-function \code{CovControlMest}.
}
\section{Slots}{
 \describe{
    \item{\code{r}:}{a numeric value specifying the required 
    breakdown point.  Allowed values are between 
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45}
    \item{\code{arp}:}{a numeric value specifying the asympthotic 
    rejection point, i.e. the fraction of points receiving zero 
    weight (see Rocke (1996)).  Default is \code{0.05} }
    \item{\code{eps}:}{a numeric value specifying the 
    relative precision of the solution of the M-estimate. 
    Defaults to \code{1e-3}}
    \item{\code{maxiter}:}{maximum number of iterations allowed 
    in the computation of the M-estimate. Defaults to 120 }
  }
}
\section{Extends}{
Class \code{"CovControl"}, directly.
}
\section{Methods}{
\describe{
    \item{estimate}{\code{signature(obj = "CovControlMest")}: the generic
    function \code{estimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovMest} passing it the control
    object and will return the obtained \code{CovRobust} object}
}

}
%\references{ ~put references to the literature/web site here ~ }
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{ ~~further notes~~ }

\seealso{
}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMest", r=0.4)
    ctrl2 <- CovControlMest(r=0.4)

    data(hbk)
    CovMest(hbk, control=ctrl1)

}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
