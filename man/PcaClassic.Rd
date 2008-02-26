\name{PcaClassic}
\alias{PcaClassic}
\alias{PcaClassic.formula}
\alias{PcaClassic.default}
\title{ Principal Components Analysis }
\description{
  Performs a principal components analysis and returns the results as an object of class PcaClassic (aka constructor).
}
\usage{
PcaClassic(x, ...)
\method{PcaClassic}{default}(x, k = 0, kmax = 10, trace=FALSE, \dots)
\method{PcaClassic}{formula}(formula, data = NULL, subset, na.action, \dots)
}
\arguments{
  \item{formula}{a formula with no response variable, referring only to
    numeric variables.}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{\dots}{arguments passed to or from other methods.}
  \item{x}{a numeric matrix (or data frame) which provides
    the data for the principal components analysis.}
  \item{k}{number of principal components to compute. If \code{k} is missing, 
    or \code{k = 0}, the algorithm itself will determine the number of 
    components by finding such \code{k} that \eqn{l_k/l_1 >= 10.E-3} and 
    \eqn{\Sigma_{j=1}^k l_j/\Sigma_{j=1}^r l_j >= 0.8}. 
    It is preferable to investigate the scree plot in order to choose the number 
    of components and the run again. Default is \code{k=0}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
%\details{}
\value{
  An S4 object of class \code{\link{PcaClassic-class}} which is a subclass of the 
  virtual class \code{\link{Pca-class}}.    
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\note{ This function can be seen as a wrapper arround prcomp() from {stats} which
returns the results of the PCA in a class compatible with the object model for robust PCA. 
}
\seealso{
     \code{\link{Pca-class}}, \code{\link{PcaClassic-class}}, 
}
%\examples{}
\keyword{robust}
\keyword{multivariate}
