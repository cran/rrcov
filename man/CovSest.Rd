\encoding{latin1}
\name{CovSest}
\alias{CovSest}
\title{ S Estimates of Multivariate Location and Scatter }
\concept{High breakdown point}
\description{
    Computes S-Estimates of multivariate location and scatter based on Tukey's
    biweight function using a fast algorithm similar to the one proposed by
    Salibian-Barrera and Yohai (2006) for the case of regression. Alternativley, the
    Ruppert's SURREAL algorithm, bisquare or Rocke type estimation can be used.
}
\usage{
    CovSest(x, bdp = 0.5, arp = 0.1, eps = 1e-5, maxiter = 120,
        nsamp = 500, seed = NULL, trace = FALSE, tolSolve = 1e-14,
        method = c("sfast", "surreal", "bisquare", "rocke", "suser"), control,
        t0, S0, initcontrol)
}
\arguments{
  \item{x}{ a matrix or data frame. }
  \item{bdp}{a numeric value specifying the required
    breakdown point.  Allowed values are between
    \code{(n - p)/(2 * n)} and 1 and the default is \code{bdp=0.5}.
  }
  \item{arp}{a numeric value specifying the asympthotic
    rejection point (for the Rocke type S estimates),
    i.e. the fraction of points receiving zero
    weight (see Rocke (1996)).  Default is \code{arp=0.1}.
  }
  \item{eps}{a numeric value specifying the
    relative precision of the solution of the S-estimate
    (bisquare and Rocke type). Default is to \code{eps=1e-5}.
  }
  \item{maxiter}{maximum number of iterations allowed
    in the computation of the S-estimate (bisquare and Rocke type).
    Default is \code{maxiter=120}.
  }
  \item{nsamp}{the number of random subsets considered. The default is different for the different methods: 
  (i) for \code{sfast} it is \code{nsamp = 20}, 
  (ii) for \code{surreal} it is \code{nsamp = 600*p} and 
  (iii) for \code{bisquare} or  \code{rocke} it is \code{nsamp = 500}.}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.}
  \item{method}{ Which algorithm to use: 'sfast'=C implementation of FAST-S, 'surreal'=SURREAL, 
  'bisquare', 'rocke'. The method 'suser' currently calls the R implementation of FAST-S 
  but in the future will allow the user to supply own \code{rho} function.}
  \item{control}{ a control object (S4) of class \code{\link{CovControlSest-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
  \item{t0}{ optional initial HBDP estimate for the center }
  \item{S0}{ optional initial HBDP estimate for the covariance matrix }
  \item{initcontrol}{ optional control object to be used for computing the initial HBDP estimates }
}
\details{
    Computes multivariate S-estimator of location and scatter. The computation will be
    performed by one of the following algorithms:
    \describe{
        \item{FAST-S}{An algorithm similar to the one proposed by Salibian-Barrera and Yohai (2006) for the case of regression}
        \item{SURREAL}{Ruppert's SURREAL algorithm when \code{method} is set to 'surreal'}
        \item{BISQUARE}{Bisquare S-Estimate with \code{method} set to 'bisquare'}
        \item{ROCKE}{Rocke type S-Estimate with \code{method} set to 'rocke'}
    }
    Except for the last algorithm, \emph{ROCKE}, all other use Tukey biweight loss function. 
    The tuning parameters used in the loss function (as determined by bdp) are 
    returned in the slots \code{cc} and \code{kp} of the result object. They can be computed 
    by the internal function \code{.csolve.bw.S(bdp, p)}.
}
\value{
  An S4 object of class \code{\link{CovSest-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  H.P. Lopuhaä (1989) On the Relation between S-estimators and M-estimators of
  Multivariate Location and Covariance.
  \emph{Annals of Statistics} \bold{17} 1662--1683.

  D. Ruppert (1992) Computing S Estimators for Regression and Multivariate
  Location/Dispersion.
  \emph{Journal of Computational and Graphical Statistics} \bold{1} 253--270.

  M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression
  estimates, \emph{Journal of Computational and Graphical Statistics}, \bold{15},
  414--427.

  R. A. Maronna, D. Martin and V. Yohai (2006). \emph{Robust Statistics: Theory and Methods}.
  Wiley, New York.

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  URL \url{http://www.jstatsoft.org/v32/i03/}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at},
Matias Salibian-Barrera \email{matias@stat.ubc.ca} and
Victor Yohai \email{vyohai@dm.uba.ar}.  See also the code from
Kristel Joossens, K.U. Leuven, Belgium and Ella Roelant, Ghent University, Belgium.
}
%\note{}
%\seealso{}
\examples{

library(rrcov)
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovSest(hbk.x)

## the following four statements are equivalent
c0 <- CovSest(hbk.x)
c1 <- CovSest(hbk.x, bdp = 0.25)
c2 <- CovSest(hbk.x, control = CovControlSest(bdp = 0.25))
c3 <- CovSest(hbk.x, control = new("CovControlSest", bdp = 0.25))

## direct specification overrides control one:
c4 <- CovSest(hbk.x, bdp = 0.40,
             control = CovControlSest(bdp = 0.25))
c1
summary(c1)
plot(c1)

## Use the SURREAL algorithm of Ruppert
cr <- CovSest(hbk.x, method="surreal")
cr

## Use Bisquare estimation
cr <- CovSest(hbk.x, method="bisquare")
cr

## Use Rocke type estimation
cr <- CovSest(hbk.x, method="rocke")
cr


}
\keyword{robust}
\keyword{multivariate}
