\name{rrcov.control}
\alias{rrcov.control}
\alias{rrcov.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control object for the estimation parameters }
\description{
     Useful for passing the estimation options as parameters to the estimation functions
}
\usage{
rrcov.control(alpha=1/2, nsamp=500, seed=0, print.it=FALSE, adjust=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{ This parameter controls the size of the subsets over which the determinant is
minimized, i.e. \code{alpha*n} observations are used for computing the determinant. Allowed values are between 0.5 and 1 and the default is 0.5. }
  \item{nsamp}{ number of subsets used for initial estimates. Default is \code{nsamp = 500}}
  \item{seed}{ starting value for random generator. Default is \code{seed = 0}}
  \item{print.it}{ whether to print intermediate results. Default is \code{print.it = FALSE}}
  \item{adjust}{ whether to perform intercept adjustment at each step. This could be quite 
time consuming, therefore the default is \code{adjust = FALSE} }
}
\details{
For details about the estimation options see the corresponding estimation functions.
}
\value{
  A list with components, as the parameters passed  by the invocation

}
\references{ 
}

% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{}
\examples{

data(hbk)
data(brain)
ctrl <- rrcov.control(alpha=0.75, print.it=TRUE)
covMcd(hbk.x, control = ctrl)
covMcd(log(brain), control = ctrl)

}
\keyword{robust}
\keyword{multivariate}
