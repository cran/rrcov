\name{SummaryCov-class}
\docType{class}
\alias{SummaryCov-class}
\alias{getCenter,SummaryCov-method}
\alias{getCov,SummaryCov-method}
\alias{getDistance,SummaryCov-method}
\alias{getEvals,SummaryCov-method}
\alias{isClassic,SummaryCov-method}
\alias{show,SummaryCov-method}

\title{Class "SummaryCov" - summary of "Cov" objects  }

\description{The "Cov" object plus some additional summary information}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummaryCov", ...)},
but most often by invoking 'summary' on a "Cov" object. They contain values 
meant for printing by 'show'.
}
\section{Slots}{
    \describe{
        \item{\code{covobj}:}{Object of class \code{"Cov"}}
        \item{\code{evals}:}{eigenvalues of the covariance or correlation matrix }
  }
}
\section{Methods}{
  \describe{
    \item{getCenter}{\code{signature(obj = "SummaryCov")}: location vector }
    \item{getCov}{\code{signature(obj = "SummaryCov")}: covariance matrix }
    \item{getDistance}{\code{signature(obj = "SummaryCov")}: vector of distances }
    \item{getEvals}{\code{signature(obj = "SummaryCov")}: vector of eignevalues }
    \item{isClassic}{\code{signature(obj = "SummaryCov")}: is the estimate a classic one }
    \item{show}{\code{signature(object = "SummaryCov")}: display the object }
     }
}
%\references{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}
\seealso{
  \code{\link{Cov-class}}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
cv <- Cov(hbk.x)
cv
summary(cv)
}
\keyword{classes}