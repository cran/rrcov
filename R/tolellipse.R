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
#
#   I would like to thank Peter Filtzmoser for providing the initial code of 
#   this function.


tolellipse <- function(x, 
                        mcd, 
                        cutoff, 
                        id.n,
                        classic=FALSE,
                        tol.inv=1e-07) {

##@bdescr
## Tolerance Ellipse Plot: 
##    Plots the 97.5% tolerance ellipse of the bivariate data set (x).
##    The ellipse is defined by those data points whose distance (dist) 
##    is equal to the squareroot of the 97.5% chisquare quantile with 
##    2 degrees of freedom. 

##@edescr
##
##@in  x                 : [matrix] A data.frame or matrix, n > 2*p
##@in  mcd               : [mcd object] An object of type mcd - its attributes 
##                                      center and cov will be used
##@in  cutoff            : [number] Distance needed to flag data points outside the ellipse 
##@in  outflag           : [logical] Whether to print the labels of the outliers 
##@in  tol.inv           : [number] tolerance to be used for computing the inverse see 'solve'.
##                                  defaults to 1e-7 

    
    ellips <- function(x, y, loc, cov) {
        # calculates a 97,5% ellipsoid 
        # input: data set, location and covariance estimate, cutoff
    
        dist <- sqrt(qchisq(0.975, 2))
        A <- solve(cov)
        lambda1 <- max(eigen(A)$values)
        lambda2 <- min(eigen(A)$values)
        eigvect <- eigen(A)$vectors[, order(eigen(A)$values)[2]]
        z <- seq(0, 2 * pi, 0.01)
        z1 <- dist/sqrt(lambda1) * cos(z)
        z2 <- dist/sqrt(lambda2) * sin(z)
        alfa <- atan(eigvect[2]/eigvect[1])
        r <- matrix(c(cos(alfa),  - sin(alfa), sin(alfa), cos(alfa)), ncol = 2)
        z <- t(t(cbind(z1, z2) %*% r) + loc)    #   xmin <- min(x, z[, 1])
    
    #   xmax <- max(x, z[, 1])
    #   ymin <- min(y, z[, 2])
    #   ymax <- max(y, z[, 2])
    #   print(xmin)
    #   print(xmax)
    #   print(ymin)
    #   print(ymax)
    #   plot(x, y, xlim = c(xmin, xmax), ylim = c(ymin, ymax), type = "n")
    #   points(z[, 1], z[, 2], type = "l")
    #   points(x, y)
        z
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
        cutoff <- sqrt(qchisq(0.975, 2))

    if(p != 2)
        stop("Dimension must be 2!")
    
    if(missing(mcd))
        mcd <- covMcd(x)

    if(length(mcd$center)  == 0 ||  length(mcd$cov) == 0)
        stop(message = "Invalid mcd object: attributes center and cov missing!")

    x.loc <- mcd$center
    x.cov <- n/(n - 1) * mcd$cov
    z1 <- ellips(x[, 1], x[, 2], loc = apply(x, 2, mean), cov = n/(n - 1) * cov.wt(x)$cov)
    z2 <- ellips(x[, 1], x[, 2], loc = x.loc, cov = x.cov)
    x1 <- c(min(x[, 1], z1[, 1], z2[, 1]), max(x[,1],z1[,1], z2[,1]))
    y1 <- c(min(x[, 2], z1[, 2], z2[, 2]), max(x[,2],z1[,2], z2[,2]))
    
    md <- sqrt(mahalanobis(x,apply(x,2,mean),cov(x), tol.inv=tol.inv))
    rd <- sqrt(mahalanobis(x,mcd$center,mcd$cov, tol.inv=tol.inv))
    
    if(classic)
        par(mfrow = c(1, 2))
    else
        par(mfrow = c(1, 1))
    
    if(missing(id.n))
        id.n <- length(which(rd>cutoff))          
    ind <- sort(rd, index.return=TRUE)$ix
    ind <- ind[(n-id.n+1):n]

##  1. Robust tollerance
##  define the plot, plot a box, plot the "good" points, 
##  plot the outliers either as points or as numbers depending on outflag,
##  plot the ellipse, write a title of the plot
    plot(x[, 1], x[, 2], xlim = x1, ylim = y1, xlab = "", ylab = "", type = "p")
    box()
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    text(x[ind, 1] + xrange/50, x[ind, 2], ind)

    points(z2[, 1], z2[, 2], type = "l")
    title(main = "ROBUST TOLERANCE \n    ELLIPSE (97.5%)")

##  2. Classical tollerance
    if(classic){
        plot(x[, 1], x[, 2], xlim = x1, ylim = y1, xlab = "", ylab = "", type = "p")
        box()
        
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(x[ind, 1] + xrange/50, x[ind, 2], ind)
        
        points(z1[, 1], z1[, 2], type = "l")
        title(main = "CLASSICAL TOLERANCE \n    ELLIPSE (97.5%)")
    }
    invisible()
}
