##  rrcov : Scalable Robust Estimators with High Breakdown Point
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
##
##  I would like to thank Peter Filtzmoser for providing the initial code of 
##  some of these functions.
##

covPlot <- function(x, 
                   plotname = c("dd","distance","qqchi2","tolellipse"),
                   classic = FALSE,
                   mcd, 
                   cutoff, 
                   nid,
                   tol.inv = 1e-7){

##@bdescr
##  Make plots based on the covariance structure of a data set:
##  dd       -  distance-distance plot: Robust distances versus 
##              Mahalanobis distances
##  distance -  a plot of the robust distances
##  qqchi2   -  a qq-plot of the robust distances versus the 
##              quantiles of the chi-squared distribution
##  tolellipse- a tolerance ellipse
##
## Distance Plot: 
## Draw a Distance-Distance Plot: Plots the robust distances
## versus the classical Mahalanobis distances as introduced by
## Rousseeuw, P. J., and van Zomeren, B. C. (1990). Unmasking
## Multivariate Outliers and Leverage Points. Journal of the American
## Statistical Association, 85, 633-639.
##
## The dashed line is the set of points where the robust distance is 
## equal to the classical distance.
## The horizontal and vertical dotted lines are drawn at values equal cutoff
## which defaults to square root of the 97.5% quantile of a chi-squared 
## distribution with p degrees of freedom. Points beyond these lines can 
## be considered outliers.
## 
##@edescr
##
##@in  x                 : [matrix] A data.frame or matrix, n > 2*p
##@in  plotname          : [character] A plot option, one of:
##                            classic: index plot of the classical mahalanobis distances
##                            robust:  index plot of the robust mahalanobis distances
##                            dd:      distance-distance plot
##                            index:   parallel index plot of classical and robust distances
##                            all:     all three plots
##                          default is "all"
##@in  classic           : [logical] If true the classical plot will be displayed too
##                                   default is classic=FALSE 
##@in  mcd               : [mcd object] An object of type mcd - its attributes 
##                                      center and cov will be used
##@in  cutoff            : [number] The cutoff value for the distances 
##@in  nid               : [number] number of observations to be identified with a label.
##                                  Defaults to the number of observations with distance
##                                  larger than cutoff 
##@in  tol.inv           : [number] tolerance to be used for computing the inverse - see 'solve'.
##                                  defaults to 1e-7 

# NOTE: The default tolerance 1e-7, will not work for some example 
#       data sets, like milk or aircraft


mydistplot <- function(x, cutoff, classic = FALSE, nid){
##  Index Plot:
##  Plot the vector x (robust or mahalanobis distances) against 
##  the observation indexes. Identify by a label the nid
##  observations with largest value of x. If nid is not supplied, 
##  calculate it as the number of observations larger than cutoff.
##  Use cutoff to draw a horisontal line.
##  Use classic=FALSE/TRUE to choose the label of the vertical axes

    n <- length(x)
    if(missing(nid))
        nid <- length(which(x>cutoff))
    if(classic)
        ylab="Square Root of Mahalanobis distance"
    else
        ylab="Square Root of Robust distance"
    plot(x, ylab=ylab, xlab="Index", type="p")
    label(1:n, x, nid)
    abline(h=cutoff)
}

myddplot <- function(md, rd, cutoff, nid){
##  Distance-Distance Plot:
##  Plot the vector y=rd (robust distances) against 
##  x=md (mahalanobis distances). Identify by a label the nid
##  observations with largest rd. If nid is not supplied, calculate
##  it as the number of observations larger than cutoff. Use cutoff
##  to draw a horisontal and a vertical line. Draw also a dotted line
##  with a slope 1.
    n <- length(md)
    if(missing(nid))
        nid <- length(which(rd>cutoff))
    xlab <- "Mahalanobis distance"
    ylab <- "Robust distance"
    plot(md, rd, xlab=xlab, ylab=ylab, type="p")
    label(md,rd,nid)
    abline(0, 1, lty=2)
    abline(v=cutoff)
    abline(h=cutoff)
}

qqplot <- function(x, p, cutoff, classic=FALSE, nid){
##  Chisquare QQ-Plot:
##  Plot the vector x (robust or mahalanobis distances) against 
##  the square root of the quantiles of the chi-squared distribution
##  with p degrees of freedom.
##  Identify by a label the nid observations with largest value of x.
##  If nid is not supplied, calculate it as the number of observations
##  larger than cutoff.
##  Use classic=FALSE/TRUE to choose the label of the vertical axes


    ##  parameters and preconditions     

    n <- length(x)

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(missing(nid))
        nid <- length(which(x>cutoff))

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

    x <- sort(x, index.return=TRUE)
    ix <- x$ix
    x <- x$x

    if(classic)
        ylab="Mahalanobis distance"
    else
        ylab="Robust distance"

    plot(qq, x, xlab="Square root of the quantiles of the chi-squared distribution", ylab=ylab, type="p")
    if(nid > 0){
        ind <- (n-nid+1):n
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(qq[ind] + xrange/50, x[ind], ix[ind])
    }
    abline(0, 1, lty=2)
}

label <- function(x, y, nid=3){
    if(nid > 0) {
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-nid+1):n]
        text(x[ind] + xrange/50, y[ind], ind)
    }
}
    
    ##  parameters and preconditions     

    if(is.vector(x) || is.matrix(x)) {
        if(!is.numeric(x))
            stop(message = "x is not a numeric dataframe or matrix.")
    }else if(is.data.frame(x)) {
        if(!all(sapply(x,data.class) == "numeric"))
            stop(message = "x is not a numeric dataframe or matrix.")
    }
    
    n <- dim(x)[1]
    p <- dim(x)[2]

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(missing(mcd))
        mcd <- covMcd(x)

    if(length(mcd$center)  == 0 ||  length(mcd$cov) == 0)
        stop(message = "Invalid mcd object: attributes center and cov missing!")

    if(length(mcd$center)  != p)
        stop(message = "Data set and provided center have different dimensions!")

    md <- mahalanobis(x, apply(x,2,mean), var(x), tol.inv=tol.inv)
    md <- sqrt(md)

    plotname <- match.arg(plotname)
    rd <- mahalanobis(x, mcd$center, mcd$cov, tol.inv=tol.inv)
    rd <- sqrt(rd)

    if(!classic || plotname == "dd")    
        par(mfrow=c(1,1), pty="m")
    else
        par(mfrow=c(1,2), pty="m")

    if(plotname == "distance"){    
        mydistplot(rd, cutoff, nid=nid)                     # index plot of mahalanobis distances
        if(classic)
            mydistplot(md, cutoff, classic=TRUE, nid=nid)   # index plot of robust distances
    }

    if(plotname == "dd"){    
        myddplot(md, rd, cutoff=cutoff, nid=nid)    # distance-distance plot
    }
    if(plotname == "qqchi2"){    
        qqplot(rd, p, cutoff=cutoff, nid=nid)       # qq-plot of the robust distances versus the 
                                                    # quantiles of the chi-squared distribution
        if(classic)
            qqplot(md, p, cutoff=cutoff, classic=TRUE, nid=nid)
                                                    # qq-plot of the mahalanobis distances
    }
    if(plotname == "tolellipse"){    
       tolellipse(x, mcd=mcd, cutoff=cutoff, nid=nid, classic=classic, tol.inv=tol.inv)
                                                # qq-plot of the robust distances versus the 
                                                # quantiles of the chi-squared distribution
    }
}

ddplot <- function(x,...){
    covPlot(x, plot="dd", ...)
}

distplot <- function(x,...){
    covPlot(x, plot="distance", ...)
}

chi2qqplot <- function(x,...){
    covPlot(x, plot="qqchi2", ...)
}

ellipse <- function(x,...){
    covPlot(x, plot="tolellipse", ...)
}
