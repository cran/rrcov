\name{ltsReg}
\alias{ltsReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust regression with high breakdown point }
\description{
Carries out least trimmed squares (LTS) regression.
}
\usage{
ltsReg(x, y, intercept=TRUE, alpha=NULL, nsamp=500, adjust=FALSE, mcd=TRUE, qr.out=FALSE, yname=NULL, seed=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{ a matrix or data frame containing the explanatory variables. }
  \item{y}{ the response: a vector of length the number of rows of 'x'. }.
  \item{intercept}{ if true, a model with constant term will be itted; otherwise no constant term will be included. Default is \code{intercept = TRUE}  }
  \item{alpha}{ the percentage of squared residuals whose sum will be minimized. Its default value is 0.5. 
  In general, alpha must be a value between 0.5 and 1.
  }
  \item{nsamp}{ number of subsets used for initial estimates. Default is \code{nsamp = 500} }
  \item{adjust}{ whether to perform intercept adjustment at each step. This could be quite time consuming, therefore the default is \code{adjust = FALSE} }
  \item{mcd}{ whether to compute robust distances using Fast-MCD.}
  \item{qr.out}{ whether to return the QR decomposition. Default is \code{qr.out = FALSE}}
  \item{yname}{ the name of the dependent variable. Default is \code{yname = NULL}}
  \item{seed}{ starting value for random generator. Default is \code{seed = 0}}
}
\details{
 The LTS regression method minimizes the sum of the h smallest squared 
residuals, where h must be at least half the number of observations. The 
default value of h is roughly 0.5n where n is the total number of observations, 
but the user may choose any value between n/2 and n. The computations are performed 
using the Fast LTS algorithm proposed by Rousseeuw and Van Driessen (1999).
}
\value{
  An object of class '"lts"'. This is a list with components

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
  \item{rsquared}{
    robust version of R squared. 
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
  \item{raw.coefficients}{
    vector of raw coefficient estimates (including the intercept,when intercept=TRUE).
  }
  \item{raw.scale}{
    scale estimate of the raw residuals.
  }
  \item{raw.resid}{
    vector like y containing the raw residuals from the regression.
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
}

\examples{

data(heart)
ltsReg(heart.x, heart.y)

}
\keyword{robust}
