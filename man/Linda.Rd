\name{Linda}
\alias{Linda}
\alias{Linda.formula}
\alias{Linda.default}
\title{ Robust Linear Discriminant Analysis }
\description{
    Robust linear discriminant analysis based on MCD and returns 
    the results as an object of class \code{Linda} (aka constructor). 
}
\usage{
Linda(x, ...)

\method{Linda}{default}(x, grouping, prior = proportions, tol = 1.0e-4,
                 method = c("mcd", "mcdA", "mcdB", "mcdC", "fsa"),
                 alpha=0.5, ...)
}
\arguments{
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{method}{method}
  \item{alpha}{this parameter measures the fraction of outliers the algorithm should
    resist. In MCD alpha controls the size of the subsets over which the 
    determinant is minimized, i.e. alpha*n observations are used for 
    computing the determinant. Allowed values are between 0.5 and 1 
    and the default is 0.5.}
  \item{\dots}{arguments passed to or from other methods}
}
\details{
  details
}
\value{
  Returns an S4 object of class \code{Linda}
}
\references{ 
    Hawkins, D.M. and McLachlan, G.J. (1997) High-Breakdown Linear Discriminant Analysis, 
    \emph{Journal of the American Statistical Association}, \bold{92}, 136--143.

    Todorov V. (2007) Robust selection of variables in linear discriminant analysis, 
        \emph{Statistical Methods and Applications}, \bold{15}, 395--407, doi:10.1007/s10260-006-0032-6.

    Todorov, V. and Pires, A.M. (2007) Comparative Performance of Several Robust 
        Linear Discriminant Analysis Methods. 
        \emph{REVSTAT Statistical Journal}, \bold{5}, p 63--83.
        URL \url{www.ine.pt/revstat/pdf/rs070104.pdf}.

    Todorov V and Filzmoser P (2009),
      An Object Oriented Framework for Robust Multivariate Analysis.
      \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
      URL \url{http://www.jstatsoft.org/v32/i03/}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\section{Warning }{Still an experimental version!}
\seealso{  \code{\link[rrcov]{CovMcd}} }
\examples{
## Example anorexia
library(MASS)
data(anorexia)

## start with the classical estimates
lda <- LdaClassic(Treat~., data=anorexia)
predict(lda)@classification

## try now the robust LDA with the default method (MCD with pooled whitin cov matrix)
rlda <- Linda(Treat~., data= anorexia)
predict(rlda)@classification

## try the other methods
Linda(Treat~., data= anorexia, method="mcdA")
Linda(Treat~., data= anorexia, method="mcdB")
Linda(Treat~., data= anorexia, method="mcdC")

## try the Hawkins&McLachlan method
## use the default method
grp <- anorexia[,1]
grp <- as.factor(grp)
x <- anorexia[,2:3]
Linda(x, grp, method="fsa")

}
\keyword{robust}
\keyword{multivariate}
