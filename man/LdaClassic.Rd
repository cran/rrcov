\name{LdaClassic}
\alias{LdaClassic}
\alias{LdaClassic.formula}
\alias{LdaClassic.default}
\title{ Linear Discriminant Analysis }
\description{
  Performs a linear discriminant analysis and returns the results as an object of class \code{LdaClassic} (aka constructor).
}
\usage{
LdaClassic(x, ...)

%%\method{LdaClassic}{formula}(formula, data = NULL, subset, na.action, \dots)
\method{LdaClassic}{default}(x, grouping, prior = proportions, tol = 1.0e-4, ...)
}
\arguments{
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{\dots}{arguments passed to or from other methods.}  
}
%\details{}
\value{
  Returns an S4 object of class \code{LdaClassic}
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\seealso{
     \code{\link{Lda-class}}, \code{\link{LdaClassic-class}}, 
}
%\examples{}
\keyword{robust}
\keyword{multivariate}
