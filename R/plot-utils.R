## Internal class defining a default legend parameters to be used 
##  in plots with both robust and classical lines in the same panel
##  (e.g. scree plot, tolerance ellipses)
setClass(".Legend", representation( leg = "logical",        # whether to draw a legend
                                    pos = "character",      # position of the legend
                                    txt = "vector",         # texts
                                    col = "vector",         # colors
                                    lty = "vector",         # line styles
                                    pch = "vector"),        # characters
                    prototype = list(leg = TRUE,
                                     pos = "topright",
                                     txt = c("robust", "classical"),
                                     col = c("red", "blue"),
                                     lty = c("solid", "dotted"),
                                     pch = c(1,24)))


.myscreeplot <- function(rcov, ccov) {
    er <- NULL
    ec <- NULL
    ee <- NULL
    if(!missing(ccov))
        ee <- ec <- getEvals(ccov)
    if(!missing(rcov))
        ee <- er <- getEvals(rcov)
    if(is.null(ee))
        stop("Both parameters are NULL")
        
    leg <- new(".Legend")        
    eall <- if(!is.null(ec) && !is.null(er))
                c(er, ec)
            else if(!is.null(ec))
                ec
            else er
                
    ylim <- c(min(eall), max(eall))
    
    plot(ee, ylim=ylim, ylab="Eigenvalues", xlab="Index", type="n")
    if(!is.null(ec) && !is.null(er))
        legend(leg@pos, leg@txt, pch = leg@pch, lty = leg@lty, col = leg@col)
    
    if(!is.null(er))
        lines(er, type="o", pch = leg@pch[1], lty = leg@lty[1], col=leg@col[1])
    if(!is.null(ec))
        lines(ec, type="o", pch = leg@pch[2], lty = leg@lty[2], col=leg@col[2])
    
    title(main="Scree plot")
}

.myddplot <- function(md, rd, cutoff, id.n) {
    ##  Distance-Distance Plot:
    ##  Plot the vector y=rd (robust distances) against
    ##  x=md (mahalanobis distances). Identify by a label the id.n
    ##  observations with largest rd. If id.n is not supplied, calculate
    ##  it as the number of observations larger than cutoff. Use cutoff
    ##  to draw a horisontal and a vertical line. Draw also a dotted line
    ##  with a slope 1.
    n <- length(md)
    if(missing(id.n))
        id.n <- length(which(rd>cutoff))
    xlab <- "Mahalanobis distance"
    ylab <- "Robust distance"
    plot(md, rd, xlab=xlab, ylab=ylab, type="p")
    .label(md,rd,id.n)
    abline(0, 1, lty=2)
    abline(v=cutoff)
    abline(h=cutoff)

    title(main="Distance-Distance Plot")
}

.mydistplot <- function(x, cutoff, classic = FALSE, id.n, ylim=NULL) {
    ## VT::10.11.2007 - change "Squared Robust distance" to "Robust distance"
    ## VT::10.11.2007 - Add parameter ylim to make possible robust and 
    ##  classical plot to use the same y-scale

    ##  Index Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the observation indexes. Identify by a label the id.n
    ##  observations with largest value of x. If id.n is not supplied,
    ##  calculate it as the number of observations larger than cutoff.
    ##  Use cutoff to draw a horisontal line.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes

    n <- length(x)
    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    ylab <- paste(if(classic) "Mahalanobis" else "Robust", "distance")
    plot(x, ylab=ylab, xlab="Index", type="p", ylim=ylim)
    .label(1:n, x, id.n)
    abline(h=cutoff)

    title(main="Distance Plot")
}

.qqplot <- function(x, p, cutoff, classic=FALSE, id.n) {
    ##  Chisquare QQ-Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the square root of the quantiles of the chi-squared distribution
    ##  with p degrees of freedom.
    ##  Identify by a label the id.n observations with largest value of x.
    ##  If id.n is not supplied, calculate it as the number of observations
    ##  larger than cutoff.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes


    ##  parameters and preconditions

    n <- length(x)

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

    x <- sort(x, index.return=TRUE)
    ix <- x$ix
    x <- x$x

    if(classic)
        ylab="Mahalanobis distance"
    else
        ylab="Robust distance"
    
    ## xlab <- "Square root of the quantiles of the chi-squared distribution"
    xlab <- eval(substitute(expression(paste("Sqrt of the quantiles of the ", chi^2, " distribution"))))
    plot(qq, x, xlab=xlab, ylab=ylab, type="p")
    if(id.n > 0) {
        ind <- (n-id.n+1):n
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(qq[ind] + xrange/50, x[ind], ix[ind])
    }
    abline(0, 1, lty=2)
    title(main=eval(substitute(expression(paste(chi^2, " QQ-Plot")))))
}

.label <- function(x, y, id.n=3) {
    if(id.n > 0) {
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-id.n+1):n]
        text(x[ind] + xrange/50, y[ind], ind)
    }
}

.tolellipse <-
    function(rcov, ccov, cutoff = NULL, id.n = NULL, tol = 1e-07, ...)
{

## MM: This is nothing else but a version  cluster::ellipsoidPoints() !! -- FIXME
    ellips <- function(loc, cov) {
        ## calculates a 97,5% ellipsoid
        ## input: data set, location and covariance estimate, cutoff
    
        dist <- sqrt(qchisq(0.975, 2))
        A <- solve(cov)
        eA <- eigen(A)
        ev <- eA$values
        lambda1 <- max(ev)
        lambda2 <- min(ev)
        eigvect <- eA$vectors[, order(ev)[2]]
        z <- seq(0, 2 * pi, by = 0.01)
        z1 <- dist/sqrt(lambda1) * cos(z)
        z2 <- dist/sqrt(lambda2) * sin(z)
        alfa <- atan(eigvect[2]/eigvect[1])
        r <- matrix(c(cos(alfa),  - sin(alfa), sin(alfa), cos(alfa)), ncol = 2)
        z <- t(t(cbind(z1, z2) %*% r) + loc)    #   xmin <- min(x, z[, 1])
    
        z
    }

    xlab <- ""
    ylab <- ""
    main <- "Tolerance ellipse (97.5%)"
    leg <- new(".Legend")
    
    if(missing(rcov) || is.null(rcov)){
        if(missing(ccov) || is.null(ccov))
            stop("Location and scatter matrix must be provided!")

        ## there is no robust location/scatter object
        rcov <- ccov
        leg@leg <- FALSE
    }
    if(is.null(data <- getData(rcov)))
        stop("No data provided!")
    n <- dim(data)[1]
    p <- dim(data)[2]
    if(p != 2)
        stop("Dimension must be 2!")

    r.cov <- getCov(rcov)
    r.loc <- getCenter(rcov)
    if(length(r.loc) == 0 ||  length(r.cov) == 0)
        stop("Invalid 'rcov' object: attributes center and cov missing!")
    z1 <- ellips(loc = r.loc, cov = r.cov)
    rd <- sqrt(getDistance(rcov))
    x1 <- c(min(data[, 1], z1[, 1]), max(data[,1],z1[,1]))
    y1 <- c(min(data[, 2], z1[, 2]), max(data[,2],z1[,2]))
    classic <- FALSE
            
    if(!missing(ccov) && !is.null(ccov)){
        c.cov <- getCov(ccov)
        c.loc <- getCenter(ccov)
        if(length(c.loc) == 0 ||  length(c.cov) == 0)
            stop("Invalid 'ccov' object: attributes center and cov missing!")
        classic <- TRUE
        z2 <- ellips(loc = c.loc, cov = c.cov)
        md <- sqrt(getDistance(ccov))
        x1 <- c(min(data[, 1], z1[, 1], z2[, 1]), max(data[,1],z1[,1], z2[,1]))
        y1 <- c(min(data[, 2], z1[, 2], z2[, 2]), max(data[,2],z1[,2], z2[,2]))
    }

    ## Note: the *calling* function may pass a 'missing' value
    if(missing(cutoff) || is.null(cutoff))
        cutoff <- sqrt(qchisq(0.975, df = 2))
    if(missing(id.n) || is.null(id.n))
        id.n <- sum(rd > cutoff)

    ind <- sort(rd, index.return=TRUE)$ix
    ind <- ind[(n-id.n+1):n]

##  1. Robust tolerance
##  define the plot, plot a box, plot the "good" points,
##  plot the outliers either as points or as numbers depending on outflag,
##  plot the ellipse, write a title of the plot
    plot(data, xlim = x1, ylim = y1, xlab = xlab, ylab = ylab)
    box()

    ## VT:::03.08.2008
    if(id.n > 0){
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(data[ind, 1] + xrange/50, data[ind, 2], ind)
    }

    if(leg@leg)
        points(z1, type = "l", lty=leg@lty[1], col=leg@col[1])
    title(main)

##  2. Classical tolerance ellipse and legend
    if(classic){
        points(z2, type = "l", lty=leg@lty[2], col=leg@col[2])
        if(leg@leg)
            legend("bottomright", leg@txt, pch=leg@pch, lty=leg@lty, col=leg@col)
    }

    invisible()
}
