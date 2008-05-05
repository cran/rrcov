\name{PcaCov}
\alias{PcaCov}
\alias{PcaCov.formula}
\alias{PcaCov.default}
\title{ Robust PCA based on a robust covariance matrix }
\description{
Robust PCA are obtained by replacing the classical covariance matrix 
by a robust covariance estimator. This can be one of the available 
in \code{rrcov} estimators, i.e. MCD, OGK, M or S estimator.
}
\usage{
PcaCov(x, ...)
\method{PcaCov}{default}(x, k = 0, kmax = ncol(x), corr=FALSE, cov.control=CovControlMcd(), 
    na.action = na.fail, trace=FALSE, \dots)
\method{PcaCov}{formula}(formula, data = NULL, subset, na.action, \dots)
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
    of components and then run again. Default is \code{k=0}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{corr}{a logical value indicating whether the calculation should use 
  the correlation matrix or the covariance matrix (the correlation matrix 
  can only be used if there are no constant variables). Default is \code{corr=FALSE}.}
  \item{cov.control}{specifies which covariance estimator to use by providing 
  a \code{\link{CovControl-class}} object. 
  The default is \code{\link{CovControlMcd-class}} which will indirectly call \code{\link{CovMcd}}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
\details{
\code{PcaCov}, serving as a constructor for objects of class \code{\link{PcaCov-class}} 
is a generic function with "formula" and "default" methods. For details see the relevant references.
}
\value{
  An S4 object of class \code{\link{PcaCov-class}} which is a subclass of the 
  virtual class \code{\link{PcaRobust-class}}. 
}
%\references{}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
%\note{}
\examples{
## PCA of the Hawkins Bradu Kass's Artificial Data
##  using all 4 variables
    data(hbk)
    pca <- PcaCov(hbk)
    pca

## Compare with the classical PCA
    prcomp(hbk)

## or  
    PcaClassic(hbk)
    
## If you want to print the scores too, use
    print(pca, print.x=TRUE)

## Using the formula interface
    PcaCov(~., data=hbk)

## To plot the results:

    plot(pca)                    # distance plot
    pca2 <- PcaCov(hbk, k=2)  
    plot(pca2)                   # PCA diagnostic plot (or outlier map)
    
## Use the standard plots available for for prcomp and princomp
    screeplot(pca)    
    biplot(pca)    
}
\keyword{robust}
\keyword{multivariate}
