\name{CovControlMcd-class}
\docType{class}
\alias{CovControlMcd-class}
\alias{estimate,CovControlMcd-method}

\title{Class 'CovControlMcd' - contains control parameters for CovMcd }
\description{ 
    This class extends the \code{CovControl} class 
    and contains the control parameters for \code{"CovMcd"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovControlMcd", ...)} 
or by calling the constructor-function \code{CovControlMcd}.
}
\section{Slots}{
     \describe{
  \item{\code{alpha}:}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{\code{nsamp}:}{ number of subsets used for initial estimates or \code{"best"}
    or \code{"exact"}.  Default is \code{nsamp = 500}.  For
    \code{nsamp="best"} exhaustive enumeration is done, as long as the
    number of trials does not exceed 5000. For \code{"exact"},
    exhaustive enumeration will be attempted however many samples are
    needed.  In this case a warning message will be displayed saying
    that the computation can take a very long time.}
  \item{\code{seed}:}{starting value for random generator. Default is \code{seed = 0}}
  %% MM (FIXME): change this to 'trace'  (as in other S functions):
  \item{\code{trace}:}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{\code{use.correction}:}{ whether to use finite sample correction factors. 
    Default is \code{use.correction=TRUE}}
  }
}
\section{Extends}{
Class \code{"CovControl"}, directly.
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(obj = "CovControlMcd")}: the generic
    function \code{estimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovMcd} passing it the control
    object and will return the obtained \code{CovRobust} object}
}}
%\references{ }
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}
%\seealso{}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMcd", alpha=0.75)
    ctrl2 <- CovControlMcd(alpha=0.75)

    data(hbk)
    CovMcd(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
