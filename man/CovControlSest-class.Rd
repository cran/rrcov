\name{CovControlSest-class}
\docType{class}
\alias{CovControlSest-class}
\alias{estimate,CovControlSest-method}
\title{Class 'CovControlSest' - contains control parameters for "CovSest"}
\description{ 
    This class extends the \code{CovControl} class 
    and contains the control parameters for CovSest
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovControlSest", ...)} 
    or by calling the constructor-function \code{CovControlSest}.
}
\section{Slots}{
     \describe{
        \item{bdp}{ a numeric value specifying the required 
            breakdown point.  Allowed values are between 
            \code{(n - p)/(2 * n)} and 1 and the default is 0.45 
        }
        \item{nsamp}{the number of random subsets considered. Default is \code{nsamp = 500}.}
        \item{\code{seed}:}{starting value for random generator. Default is \code{seed = NULL}.}
        \item{trace}{whether to print intermediate results. Default 
            is \code{trace = FALSE}.
        }  
        \item{tolSolve}{numeric tolerance to be used for inversion
            (\code{\link{solve}}) of the covariance matrix in
            \code{\link{mahalanobis}}.
        }
        \item{algo}{ Which algorithm to use: 'sfast'=FAST-S or 'surreal'=SURREAL }
     
    }
}
\section{Extends}{
    Class \code{"\linkS4class{CovControl}"}, directly.
}
\section{Methods}{
\describe{
    \item{estimate}{\code{signature(obj = "CovControlSest")}: the generic
    function \code{estimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovSest} passing it the control
    object and will return the obtained \code{CovRobust} object}
}
}
%\references{ }
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}
%\seealso{}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlSest", bdp=0.4)
    ctrl2 <- CovControlSest(bdp=0.4)

    data(hbk)
    CovSest(hbk, control=ctrl1)

}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
