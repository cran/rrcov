\name{tolellipse}
\alias{tolellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tolerance Ellipse Plot }
\description{
 Plots the 0.975 tolerance ellipse of the bivariate data set \code{x}.
The ellipse is defined by those data points whose distance 
is equal to the squareroot of the 0.975 chisquare quantile 
with 2 degrees of freedom. 
}
\usage{
tolellipse(x, mcd, cutoff, nid, classic = FALSE, tol.inv=1e-7)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a two dimensional matrix or data frame. }
  \item{mcd}{ An object of type mcd - its attributes center and cov will be used. If missing the MCD will be computed. }
  \item{cutoff}{ Distance needed to flag data points outside the ellipse  }
  \item{nid}{ Number of observations to identify by a label. If not supplied, the number of observations with distance larger than \code{cutoff} is used.  }
  \item{classic}{ Whether to plot the classical distances too. Default is \code{classic}=FALSE. }.
  \item{tol.inv}{ Tolerance to be used for computing the inverse see 'solve'. Default is \code{tol.inv = 10e-7}}
}

%\details{}
%\value{}
%\references{ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tolellipse}} 
}

\examples{

data(hbk)
mcd <- covMcd(hbk.x)       # compute mcd in advance
tolellipse(hbk.x[,1:2])    # must be a 2-dimensional data set: take the first two columns
}
\keyword{robust}
\keyword{multivariate}
