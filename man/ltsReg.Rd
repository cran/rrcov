\name{ltsReg}
\alias{ltsReg}
\alias{ltsReg.default}
\alias{ltsReg.formula}
\alias{print.lts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust regression with high breakdown point }
\description{
Carries out least trimmed squares (LTS) regression.
}
\usage{
ltsReg(x, \dots)

\method{ltsReg}{formula}(formula, data, \dots, model = TRUE, x.ret = FALSE, y.ret = FALSE)

\method{ltsReg}{default}(x, y, intercept=TRUE, alpha=NULL, nsamp=500, adjust=FALSE, mcd=TRUE, qr.out=FALSE, yname=NULL, seed=0, use.correction=TRUE, control, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{formula}{ a formula of the form 'y ~ x1 + x2 + ...'. }
  \item{data}{ data frame from which variables specified in 'formula' are to be taken. }
  \item{model, x.ret, y.ret}{ logical. If 'TRUE' the model frame, the model matrix and the response are returned, respectively. }
  \item{x}{ a matrix or data frame containing the explanatory variables. }
  \item{y}{ the response: a vector of length the number of rows of 'x'. }.
  \item{intercept}{ if true, a model with constant term will be estimated; otherwise no constant term will be included. Default is \code{intercept = TRUE}  }
  \item{alpha}{ 
    the percentage of squared residuals whose sum will be minimized. Its default value is 0.5. 
    In general, alpha must be a value between 0.5 and 1.
  }
  \item{nsamp}{ 
    number of subsets used for initial estimates or \code{"best"} or \code{"exact"}. 
    Default is \code{nsamp = 500}. 
    If \code{nsamp="best"} exhaustive enumeration is done, as far as the number of trials do not exceed 5000.
    If \code{nsamp="exact"} exhaustive enumeration will be attempted however many samples are needed. 
    In this case a warning message will be displayed saying that the computation can take a very long time. }
  \item{adjust}{ whether to perform intercept adjustment at each step. This could be quite time consuming, therefore the default is \code{adjust = FALSE} }
  \item{mcd}{ whether to compute robust distances using Fast-MCD.}
  \item{qr.out}{ whether to return the QR decomposition. Default is \code{qr.out = FALSE}}
  \item{yname}{ the name of the dependent variable. Default is \code{yname = NULL}}
  \item{seed}{ starting value for random generator. Default is \code{seed = 0}}
  \item{use.correction}{ whether to use finite sample correction factors. Default is \code{use.correction=TRUE}}
  \item{control}{ a list with estimation options - same as these provided in the 
  fucntion specification. If the control object is supplied, the parameters from it 
  will be used. If parameters are passed also in the invocation statement, they will 
  override the corresponding elements of the control object. }
  \item{\dots}{ arguments passed to or from other methods. } 
}
\details{
 The LTS regression method minimizes the sum of the h smallest squared 
residuals, where h must be at least half the number of observations. The 
default value of h is roughly 0.5n where n is the total number of observations, 
but the user may choose any value between n/2 and n. 
The LTS estimate of the error scale is given by the minimum of the objective 
function multiplied by a consistency factor and a finite sample correction 
factor - see Pison et.al. (2002) for details. The rescaling factors for the raw and final estimates are 
returned also in the vectors \code{raw.cnp2} and \code{cnp2} of length 2 respectively.
The finite sample corrections can be suppressed by setting \code{use.correction=FALSE}.
The computations are performed 
using the Fast LTS algorithm proposed by Rousseeuw and Van Driessen (1999).

The formula interface has an implied intercept term.  This can be removed by using
either \code{y ~ x - 1} or \code{y ~ 0 + x}.  See \code{\link{formula}} for more details.
}
\value{
The function \code{ltsReg} returns an object of class \code{"lts"}. 
The function \code{summary} is used to obtain and print a summary table of the results.  
The generic accessor functions \code{coefficients}, \code{fitted.values} and \code{residuals}
extract various useful features of the value returned by \code{ltsReg}.

An object of class \code{lts} is a list containing at least the following components:

  \item{crit}{
    the value of the objective function of the LTS regression method, i.e. the sum of the h smallest squared raw residuals.
  }
  \item{coefficients}{
    vector of coefficient estimates (including the intercept,when intercept=TRUE), obtained after reweighting
  }
  \item{best}{
    the best subset found and used for computing the raw estimates. The size of \code{best} is equal to \code{quan}.
  }
  \item{fitted.values}{
    vector like y containing the fitted values of the response after reweighting.
  }
  \item{residuals}{
    vector like y containing the residuals from the weighted least squares regression.
  }
  \item{scale}{
    scale estimate of the reweighted residuals. 
  }
  \item{alpha}{ 
    same as the input parameter \code{alpha}.
  }  
  \item{quan}{ 
   the number h of observations that have determined the least trimmed squares estimator
  }  
  \item{intercept}{ 
    same as the input parameter \code{intercept}. 
  }  
  \item{cnp2}{
    a vector of length two containing the consistency correction factor and the 
    finite sample correction factor of the final estimate of the error scale.
  }
  \item{raw.coefficients}{
    vector of raw coefficient estimates (including the intercept,when intercept=TRUE).
  }
  \item{raw.scale}{
    scale estimate of the raw residuals.
  }
  \item{raw.resid}{
    vector like y containing the raw residuals from the regression.
  }
  \item{raw.cnp2}{
    a vector of length two containing the consistency correction factor and the 
    finite sample correction factor of the raw estimate of the error scale.
  }
  \item{lts.wt}{
    vector like y containing weights that can be used in a weighted
    least squares. These weights are 1 for points with reasonably
    small raw residuals, and 0 for points with large raw residuals.
  }
  \item{method}{
    character string naming the method (Least Trimmed Squares).
  }
  \item{X}{
    the input data as a matrix.
  }
  \item{Y}{
    the response variable as a vector.
  }
  
}
\references{ 

 p. j. Rousseeuw (1984), Least Median of Squares Regression.
 \emph{Journal of the American Statistical Association}, \bold{79}, pp. 871-881. 
 
 P. J. Rousseeuw and A. M. Leroy (1987) 
 \emph{Robust Regression and Outlier Detection.} Wiley. 
 
 P. J. Rousseeuw and K. van Driessen (1999)
 Computing LTS Regression for Large Data Sets, 
 Technical Report, University of Antwerp, submitted
 
 P. J. Rousseeuw and K. van Driessen (1999) 
 A fast algorithm for the minimum covariance determinant estimator. 
 \emph{Technometrics}, \bold{41}, 212--223.

 Pison, G., Van Aelst, S., and Willems, G. (2002), 
 Small Sample Corrections for LTS and MCD, 
 \emph{Metrika}, \bold{55}, 111-123.
 
}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{covMcd}} 

\code{\link{summary.lts}} for summaries.
  
The generic functions \code{\link{coef}}, \code{\link{residuals}}, \code{\link{fitted}}.

}

\examples{

data(heart)
ltsReg(heart.x, heart.y)

data(stackloss)
ltsReg(stack.loss ~ ., data = stackloss)

}
\keyword{robust}
