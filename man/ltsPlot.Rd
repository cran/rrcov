\name{plot.lts}
\alias{plot.lts}
\alias{ltsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust Regression Diagnostic Plots }
\description{
    Four plots (selectable by 'which') are currently provided: a plot
of the standardized residuals versus their index, a plot
of the standardized residuals versus fitted values, a Normal Q-Q plot of the standardized residuals, and
a regression diagnostic plot (standardized residuals versus robust distances of the predictor variables).

}
\usage{
    ## S3 method for class 'lts':
    plot.lts(x, which = c("all","rqq","rindex","rfit","rdiag"), classic=FALSE, ask=(which=="all" && dev.interactive()), id.n=3, ...)
    
%    ltsPlot(x, which = c("all","rqq","rindex","rfit","rdiag"), classic=FALSE, ask=(which=="all" && dev.interactive()), id.n=3, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x }{ a \code{lts} object, typically result of \code{ltsReg}. }
  \item{which }{ Which plot to show? See Details for description of the options. Default is \code{which}="all". }.
  \item{classic }{ whether to plot the classical plots too. Default is \code{classic}=FALSE. }.
  \item{ask}{ logical; if 'TRUE', the user is \emph{ask}ed before each plot, see 'par(ask=.)'. 
Default is \code{ask = which=="all" && dev.interactive()}.  }
  \item{id.n }{ Number of observations to identify by a label starting with the most extreme. Default is \code{id.n} = 3.  }
  \item{... }{ other parameters to be passed through to plotting functions. }
  
}
\details{

 This function produces several plots based on the robust and classical regression estimates. 
Which of them to select is specified by the attribute  \code{which}.
The possible options are: 

 \code{rqq} - Normal Q-Q plot of the standardized residuals;
 
 \code{rindex} - plot of the standardized residuals versus their index;
 
 \code{rfit} - plot of the standardized residuals versus fitted values;
 
 \code{tolellipse} - regression diagnostic plot.

 The regression diagnostic Plot, introduced by
Rousseeuw and van Zomeren (1990), displays the standardized residuals versus robust distances. 
Following Rousseeuw and van Zomeren (1990), the horizontal dashed lines are located at +2.5 and -2.5 
and the vertical line is located at the upper 0.975 percent point of the chi-squared distribution 
with p degrees of freedom.

}

}
%\value{}
\references{ 
      P. J. Rousseeuw and van Zomeren, B. C. (1990). Unmasking
    Multivariate Outliers and Leverage Points. \emph{Journal of the American
    Statistical Association} \bold{85}, 633-639.
    
    P. J. Rousseeuw and K. van Driessen (1999) 
    A fast algorithm for the minimum covariance determinant estimator. 
    \emph{Technometrics} \bold{41}, 212--223.
    
 }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{covPlot}} 
}

\examples{

data(hbk)
lts <- ltsReg(hbk.x,hbk.y)
plot(lts, which="rqq")

}
\keyword{robust}
\keyword{multivariate}
