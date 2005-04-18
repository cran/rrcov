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
##  I would like to thank Peter Rousseeuw and Katrien van Driessen for 
##  providing the initial code of this function.


ltsReg <- function (x, y, 
                    intercept=TRUE, 
                    alpha=NULL, 
                    nsamp=500, 
                    adjust=FALSE, 
                    mcd=TRUE, 
                    qr.out=FALSE, 
                    yname=NULL, 
                    seed=0) 
{

    quan.f <- function(alpha, n, rk) {
        quan <- floor(2*floor((n+rk+1)/2) - n + 2*(n - floor((n+rk+1)/2)) * alpha)
        return(quan)
    }

    correctiefactor.s <- function(p, intercept = intercept, n, alpha) {
        if (intercept == TRUE) {
            p <- p - 1
        }
        if (p == 0) {
            fp.500.n <- 1 - exp(0.262024211897096) * 1/n^{
                0.604756680630497
            }
            fp.875.n <- 1 - exp(-0.351584646688712) * 1/n^{
                1.01646567502486
            }
            if ((0.5 <= alpha) && (alpha <= 0.875)) {
                fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * 
                  (alpha - 0.5)
                fp.alpha.n <- sqrt(fp.alpha.n)
            }
            if ((0.875 < alpha) && (alpha < 1)) {
                fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * 
                  (alpha - 0.875)
                fp.alpha.n <- sqrt(fp.alpha.n)
            }
        }
        else {
            if (p == 1) {
                if (intercept == TRUE) {
                  fp.500.n <- 1 - exp(0.630869217886906) * 1/n^{
                    0.650789250442946
                  }
                  fp.875.n <- 1 - exp(0.565065391014791) * 1/n^{
                    1.03044199012509
                  }
                }
                else {
                  fp.500.n <- 1 - exp(-0.0181777452315321) * 
                    1/n^{
                    0.697629772271099
                  }
                  fp.875.n <- 1 - exp(-0.310122738776431) * 1/n^{
                    1.06241615923172
                  }
                }
            }
            if (p > 1) {
                if (intercept == TRUE) {
                  coefgqpkwad875 <- matrix(c(-0.458580153984614, 
                    1.12236071104403, 3, -0.267178168108996, 
                    1.1022478781154, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefgqpkwad875) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefgqpkwad875mint.q3", 
                    "coefgqpkwad875mint.q5"))
                  coefeqpkwad500 <- matrix(c(-0.746945886714663, 
                    0.56264937192689, 3, -0.535478048924724, 
                    0.543323462033445, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefeqpkwad500) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefeqpkwad500mint.q3", 
                    "coefeqpkwad500mint.q5"))
                }
                else {
                  coefgqpkwad875 <- matrix(c(-0.251778730491252, 
                    0.883966931611758, 3, -0.146660023184295, 
                    0.86292940340761, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefgqpkwad875) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefgqpkwad875zint.q3", 
                    "coefgqpkwad875zint.q5"))
                  coefeqpkwad500 <- matrix(c(-0.487338281979106, 
                    0.405511279418594, 3, -0.340762058011, 0.37972360544988, 
                    5), ncol = 2, byrow = FALSE)
                  dimnames(coefeqpkwad500) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefeqpkwad500zint.q3", 
                    "coefeqpkwad500zint.q5"))
                }
                y1.500 <- 1 + coefeqpkwad500[1, 1] * 1/p^{
                  coefeqpkwad500[2, 1]
                }
                y2.500 <- 1 + coefeqpkwad500[1, 2] * 1/p^{
                  coefeqpkwad500[2, 2]
                }
                y1.875 <- 1 + coefgqpkwad875[1, 1] * 1/p^{
                  coefgqpkwad875[2, 1]
                }
                y2.875 <- 1 + coefgqpkwad875[1, 2] * 1/p^{
                  coefgqpkwad875[2, 2]
                }
                y1.500 <- log(1 - y1.500)
                y2.500 <- log(1 - y2.500)
                y.500 <- c(y1.500, y2.500)
                A.500 <- matrix(c(1, log(1/(coefeqpkwad500[3, 
                  1] * p^2)), 1, log(1/(coefeqpkwad500[3, 2] * 
                  p^2))), ncol = 2, byrow = TRUE)
                coeffic.500 <- solve(A.500, y.500)
                y1.875 <- log(1 - y1.875)
                y2.875 <- log(1 - y2.875)
                y.875 <- c(y1.875, y2.875)
                A.875 <- matrix(c(1, log(1/(coefgqpkwad875[3, 
                  1] * p^2)), 1, log(1/(coefgqpkwad875[3, 2] * 
                  p^2))), ncol = 2, byrow = TRUE)
                coeffic.875 <- solve(A.875, y.875)
                fp.500.n <- 1 - exp(coeffic.500[1]) * 1/n^{
                  coeffic.500[2]
                }
                fp.875.n <- 1 - exp(coeffic.875[1]) * 1/n^{
                  coeffic.875[2]
                }
            }
            if ((0.5 <= alpha) && (alpha <= 0.875)) {
                fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * 
                  (alpha - 0.5)
            }
            if ((0.875 < alpha) && (alpha <= 1)) {
                fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * 
                  (alpha - 0.875)
            }
        }
        return(1/fp.alpha.n)
    }
    correctiefactor.rew.s <- function(p, intercept = intercept, n, alpha) {
        if (intercept == TRUE) {
            p <- p - 1
        }
        if (p == 0) {
            fp.500.n <- 1 - (exp(1.11098143415027) * 1)/n^{
                1.5182890270453
            }
            fp.875.n <- 1 - (exp(-0.66046776772861) * 1)/n^{
                0.88939595831888
            }
            if ((0.5 <= alpha) && (alpha <= 0.875)) {
                fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * 
                  (alpha - 0.5)
                fp.alpha.n <- sqrt(fp.alpha.n)
            }
            if ((0.875 < alpha) && (alpha < 1)) {
                fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * 
                  (alpha - 0.875)
                fp.alpha.n <- sqrt(fp.alpha.n)
            }
        }
        else {
            if (p == 1) {
                if (intercept == TRUE) {
                  fp.500.n <- 1 - (exp(1.58609654199605) * 1)/n^{
                    1.46340162526468
                  }
                  fp.875.n <- 1 - (exp(0.391653958727332) * 1)/n^{
                    1.03167487483316
                  }
                }
                else {
                  fp.500.n <- 1 - (exp(0.6329852387657) * 1)/n^{
                    1.40361879788014
                  }
                  fp.875.n <- 1 - (exp(-0.642240988645469) * 
                    1)/n^{
                    0.926325452943084
                  }
                }
            }
            if (p > 1) {
                if (intercept == TRUE) {
                  coefqpkwad875 <- matrix(c(-0.474174840843602, 
                    1.39681715704956, 3, -0.276640353112907, 
                    1.42543242287677, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefqpkwad875) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefqpkwad875mint.q3", 
                    "coefqpkwad875mint.q5"))
                  coefqpkwad500 <- matrix(c(-0.773365715932083, 
                    2.02013996406346, 3, -0.337571678986723, 
                    2.02037467454833, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefqpkwad500) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefqpkwad500mint.q3", 
                    "coefqpkwad500mint.q5"))
                }
                else {
                  coefqpkwad875 <- matrix(c(-0.267522855927958, 
                    1.17559984533974, 3, -0.161200683014406, 
                    1.21675019853961, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefqpkwad875) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefqpkwad875zint.q3", 
                    "coefqpkwad875zint.q5"))
                  coefqpkwad500 <- matrix(c(-0.417574780492848, 
                    1.83958876341367, 3, -0.175753709374146, 
                    1.8313809497999, 5), ncol = 2, byrow = FALSE)
                  dimnames(coefqpkwad500) <- list(c("alfaq", 
                    "betaq", "qwaarden"), c("coefqpkwad500zint.q3", 
                    "coefqpkwad500zint.q5"))
                }
                y1.500 <- 1 + (coefqpkwad500[1, 1] * 1)/p^{
                  coefqpkwad500[2, 1]
                }
                y2.500 <- 1 + (coefqpkwad500[1, 2] * 1)/p^{
                  coefqpkwad500[2, 2]
                }
                y1.875 <- 1 + (coefqpkwad875[1, 1] * 1)/p^{
                  coefqpkwad875[2, 1]
                }
                y2.875 <- 1 + (coefqpkwad875[1, 2] * 1)/p^{
                  coefqpkwad875[2, 2]
                }
                y1.500 <- log(1 - y1.500)
                y2.500 <- log(1 - y2.500)
                y.500 <- c(y1.500, y2.500)
                A.500 <- matrix(c(1, log(1/(coefqpkwad500[3, 
                  1] * p^2)), 1, log(1/(coefqpkwad500[3, 2] * 
                  p^2))), ncol = 2, byrow = TRUE)
                coeffic.500 <- solve(A.500, y.500)
                y1.875 <- log(1 - y1.875)
                y2.875 <- log(1 - y2.875)
                y.875 <- c(y1.875, y2.875)
                A.875 <- matrix(c(1, log(1/(coefqpkwad875[3, 
                  1] * p^2)), 1, log(1/(coefqpkwad875[3, 2] * 
                  p^2))), ncol = 2, byrow = TRUE)
                coeffic.875 <- solve(A.875, y.875)
                fp.500.n <- 1 - (exp(coeffic.500[1]) * 1)/n^{
                  coeffic.500[2]
                }
                fp.875.n <- 1 - (exp(coeffic.875[1]) * 1)/n^{
                  coeffic.875[2]
                }
            }
            if ((0.5 <= alpha) && (alpha <= 0.875)) {
                fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * 
                  (alpha - 0.5)
            }
            if ((0.875 < alpha) && (alpha <= 1)) {
                fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * 
                  (alpha - 0.875)
            }
        }
        return(1/fp.alpha.n)
    }
    
#cat("++++++ Entering ltsReg() ...\n")    
    if (is.vector(y) || (is.matrix(y) && !is.data.frame(y))) {
        if (!is.numeric(y)) 
            stop(message = "y is not a numeric dataframe or vector.")
    }
    if ((!is.matrix(y) && !is.vector(y)) || is.data.frame(y)) {
        if ((!is.data.frame(y) && !is.numeric(y)) || (!all(sapply(y, data.class) == "numeric"))) 
            stop(message = "y is not a numeric dataframe or vector.")
    }
    
    y <- as.matrix(y)
    if (dim(y)[2] != 1) 
        stop(message = "y is not onedimensional.")
    
    if (missing(x)) {
#cat("++++++ Prepare: x is missing...\n")    
        x <- rep(1, nrow(y))
        if (is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
            if (!is.numeric(x)) 
                stop(message = "x is not a numeric dataframe or matrix.")
        }
        if ((!is.matrix(x) && !is.vector(x)) || is.data.frame(x)) {
            if ((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x, 
                data.class) == "numeric"))) 
                stop(message = "x is not a numeric dataframe or matrix.")
        }
        if (!is.matrix(x)) 
            x <- array(x, c(length(x), 1), list(names(x), deparse(substitute(x))))
        x <- as.matrix(x)
        dimny <- dimnames(y)
        dimnx <- dimnames(x)
        na.x <- !is.finite(x %*% rep(1, ncol(x)))
        na.y <- !is.finite(y)
        if (nrow(na.x) != nrow(na.y)) 
            stop("Number of observations in x and y not equal")
        ok <- !(na.x | na.y)
        y <- y[ok, , drop = FALSE]
        dy <- nrow(y)
        rownames <- dimny[[1]]
        yn <- if (!is.null(yname)) yname else dimny[[2]]
        if (!length(yn)) 
            yn <- "Y"
        storage.mode(y) <- "double"
        x <- x[ok, , drop = FALSE]
        storage.mode(x) <- "double"
        dx <- dim(x)
        if (!length(dx)) 
            stop("All observations have missing values!")
        n <- dx[1]
    }else {

#cat("++++++ Prepare: x is present...\n")    
        if (is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
            if (!is.numeric(x)) 
                stop(message = "x is not a numeric dataframe or matrix.")
        }
        if ((!is.matrix(x) && !is.vector(x)) || is.data.frame(x)) {
            if ((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x, data.class) == "numeric"))) 
                stop(message = "x is not a numeric dataframe or matrix.")
        }

        #vt:: if the data is supplied as a data.frame, the following expressions results in an error
        # as workaround convert the data.frame to a matrix
        if(is.data.frame(x))
            x <- as.matrix(x)
    
        if (!is.matrix(x)) 
            x <- array(x, c(length(x), 1), list(names(x), deparse(substitute(x))))
        x <- as.matrix(x)
        dimny <- dimnames(y)
        dimnx <- dimnames(x)
        na.x <- !is.finite(x %*% rep(1, ncol(x)))
        na.y <- !is.finite(y)
        if (nrow(na.x) != nrow(na.y)) 
            stop("Number of observations in x and y not equal")
        ok <- !(na.x | na.y)
        y <- y[ok, , drop = FALSE]
        dy <- nrow(y)
        rownames <- dimny[[1]]
        yn <- if (!is.null(yname)) 
            yname
        else dimny[[2]]
        if (!length(yn)) 
            yn <- "Y"
        storage.mode(y) <- "double"
        x <- x[ok, , drop = FALSE]
        storage.mode(x) <- "double"
        dx <- dim(x)
        if (!length(dx)) 
            stop("All observations have missing values!")
        n <- dx[1]
        constantcolom <- function(x) {
            c1 <- range(x)
            c1[1] == c1[2]
        }
        if (sum(apply(x, 2, constantcolom)) > 0) 
            stop("There is at least one constant column. Remove this column and set intercept=T")
    }

#cat("++++++ Prepare: Ready.\n")    

    dn <- dimnames(x)
    xn <- dn[[2]]
    if (!length(xn)) 
        if (dx[2] > 1) 
            xn <- paste("X", 1:dx[2], sep = "")
        else xn <- "X"
    X <- x
    dimnames(X) <- list(NULL, xn)
    y <- as.vector(y)
    
    
    if (all(x == 1)) {

#cat("++++++ A - all x == 1...\n")    
        if (length(alpha)) {
            if (alpha < 1/2) 
                stop("alpha is out of range!")
            else if (alpha > 1) 
                stop("alpha is greater than 1")
            quan <- quan.f(alpha, n, dx[2])
        }
        else {
            alpha <- 1/2
            quan <- quan.f(alpha, n, dx[2])
        }
        storage.mode(quan) <- "integer"
        initcov <- matrix(0, nrow = 1, ncol = 1)
        initmean <- matrix(0, nrow = 1, ncol = 1)
        inbest <- matrix(0, nrow = quan, ncol = 1)
        weights <- matrix(0, nrow = n, ncol = 1)
        coeff <- matrix(0, nrow = 5, ncol = 1)
        adjustcov <- matrix(0, nrow = 1, ncol = 1)
        storage.mode(initcov) <- "double"
        storage.mode(initmean) <- "double"
        storage.mode(inbest) <- "integer"
        storage.mode(weights) <- "integer"
        storage.mode(coeff) <- "double"
        storage.mode(adjustcov) <- "double"
        storage.mode(y) <- "double"
        storage.mode(seed) <- "integer"
        exactfit <- 0
        p <- 1
        xbest <- NULL
        if (alpha == 1) {
            scale <- sqrt(cov.wt(x)$cov)
            center <- as.vector(mean(x))
        }else {
            sh <- .Fortran("rffastmcd", 
                    as.matrix(y), 
                    as.integer(n), 
                    as.integer(p), 
                    as.integer(quan), 
                    nsamp = 0, 
                    initcovariance = initcov, 
                    initmean = initmean, 
                    inbest = inbest, 
                    mcdestimate = 0, 
                    weights = weights, 
                    as.integer(exactfit), 
                    coeff = coeff, 
                    kount = 0, 
                    adjustcov = adjustcov,
                    seed,
                    PACKAGE="rrcov")

            y <- as.vector(y)
            center <- as.double(sh$initmean)
            qalpha <- qchisq(quan/n, 1)
            calphainvers <- pgamma(qalpha/2, 1/2 + 1)/(quan/n)
            calpha <- 1/calphainvers
            correct <- correctiefactor.s(1, intercept = intercept, 
                n, alpha)
            scale <- sqrt(as.double(sh$initcovariance)) * sqrt(calpha) * 
                correct
            xbest <- sort(as.vector(sh$inbest))
        }
        resid <- y - center
        ans <- list()
        ans$best <- xbest
        ans$coefficients <- center
        ans$alpha <- alpha
        ans$quan <- quan
        ans$raw.resid <- resid/scale
        weights <- rep(NA, n)
        if (abs(scale) < 1e-07) {
            weights <- ifelse(abs(resid) < 1e-07, 1, 0)
            ans$scale <- ans$raw.scale <- 0
            ans$crit <- 0
            ans$coefficients <- ans$raw.coefficients <- center
        }
        if (abs(scale) >= 1e-07) {
            ans$raw.scale <- scale
            ans$raw.coefficients <- center
            quantiel <- qnorm(0.9875)
            weights <- ifelse(abs(resid/scale) <= quantiel, 1, 
                0)
            reweighting <- cov.wt(y, wt = weights)
            ans$coefficients <- reweighting$center
            ans$scale <- sqrt(sum(weights)/(sum(weights) - 1) * 
                reweighting$cov)
            resid <- y - ans$coefficients
            ans$crit <- sum(sort((y - center)^2, quan)[1:quan])
            if (sum(weights) == n) {
                cdelta.rew <- 1
                correct.rew <- 1
            }
            else {
                qdelta.rew <- qchisq(sum(weights)/n, 1)
                cdeltainvers.rew <- pgamma(qdelta.rew/2, 1/2 + 
                  1)/(sum(weights)/n)
                cdelta.rew <- sqrt(1/cdeltainvers.rew)
                correct.rew <- correctiefactor.rew.s(1, intercept = intercept, 
                  n, alpha)
            }
            ans$scale <- ans$scale * cdelta.rew * correct.rew
            quantiel <- qnorm(0.9875)
            weights <- ifelse(abs(resid/ans$scale) <= quantiel, 
                1, 0)
        }
        ans$resid <- resid/ans$scale
        ans$rsquared <- 0
        ans$residuals <- rep(NA, length(na.y))
        ans$residuals[ok] <- resid
        ans$lts.wt <- rep(NA, length(na.y))
        ans$lts.wt[ok] <- weights
        ans$intercept <- intercept
        ans$method <- paste("Univariate location and scale estimation.")
        if (abs(scale) < 1e-07) 
            ans$method <- paste(ans$method, "\nMore than half of the data are equal!")
        names(ans$coefficients) <- names(ans$raw.coefficients) <- yn
        names(ans$scale) <- names(ans$raw.scale) <- yn
        names(ans$rsquared) <- yn
        names(ans$crit) <- yn
        names(ans$residuals) <- rownames
        names(ans$lts.wt) <- rownames
        ans$X <- x
        ans$Y <- y                  # VT:: 01.09.2004 - add y to the result object
        if (length(rownames)) 
            dimnames(ans$X)[[1]] <- rownames[ok]
        else {
            xx <- seq(1, length(na.x))
            dimnames(ans$X) <- list(NULL, NULL)
            dimnames(ans$X)[[1]] <- xx[ok]
        }
        class(ans) <- "lts"
        attr(ans, "call") <- sys.call()

#cat("++++++ A - all x == 1...Ready and Return.\n")    
        return(ans)
    }

    ans <- list()
    if (intercept) {
        dx <- dx + c(0, 1)
        xn <- c(xn, "Intercept")
        x <- array(c(x, rep(1, n)), dx, dimnames = list(dn[[1]], xn))
    }
    p <- dx[2]
    if (n <= 2 * p) 
        stop("Need more than twice as many observations as variables.")
    if (length(alpha)) {
        if (alpha > 1) 
            stop("alpha is greater than 1")
        if (alpha == 1) {                                   # alpha == 1 -----------------------

#cat("++++++ B - alpha == 1...\n")    
            z <- lsfit(x, y, intercept = FALSE)
            
            # VT:: 26.12.2004 
            # Reorder the coeficients,so that the intercept moves to the beginning of the array
            # Skip this if p == 1 (i.e. p=1 and intercept=FALSE).
            # Do the same for the names and for ans$coef - see below
            if(p > 1)
                ans$raw.coefficients[2:p] <- z$coef[1:(p - 1)]
            ans$raw.coefficients[1] <- z$coef[p]
            
            ans$alpha <- alpha
            ans$quan <- quan <- n           # VT:: 01.09.2004 - bug in alpha=1 
                                            # (ans$quan was not set)
            # VT:: 26.12.2004 
            if(p > 1)
                names(ans$raw.coefficients)[2:p] <- xn[1:(p - 1)]
            names(ans$raw.coefficients)[1] <- xn[p]
            
            s0 <- sqrt((1/(n - p)) * sum(z$residuals^2))
            weights <- rep(NA, n)
#cat("++++++ B - alpha == 1... - s0=",s0,"\n")    
            if(abs(s0) < 1e-07) {
                fitted <- x %*% z$coef
                weights <- ifelse(abs(z$residuals) <= 1e-07, 1, 0)
                ans$scale <- ans$raw.scale <- 0
                ans$coefficients <- ans$raw.coefficients
            }
            else {
                ans$raw.scale <- s0
                ans$raw.resid <- ans$residuals/ans$raw.scale
                weights <- ifelse(abs(z$residuals/s0) <= qnorm(0.9875), 1, 0)
                
                # vt:: weights has to be a vector instead of a matrix - 
                #      to avoid "Error in x * wtmult : non-conformable arrays"
                # 
                weights <- as.vector(weights)
                z <- lsfit(x, y, wt = weights, intercept = FALSE)

                # VT:: 26.12.2004 
                if(p > 1)
                    ans$coefficients[2:p] <- z$coef[1:(p - 1)]
                ans$coefficients[1] <- z$coef[p]
                
                fitted <- x %*% z$coef
                ans$scale <- sqrt(sum(weights * z$residuals^2)/(sum(weights) - 
                  1))
                if (sum(weights) == n) {
                  cdelta.rew <- 1
                }
                else {
                  cdelta.rew <- (1/sqrt(1 - ((2 * n)/(sum(weights) * 
                    (1/qnorm((sum(weights) + n)/(2 * n))))) * 
                    dnorm(1/(1/(qnorm((sum(weights) + n)/(2 * 
                      n)))))))
                }
                ans$scale <- ans$scale * cdelta.rew
                weights <- ifelse(abs(z$residuals/ans$scale) <= 
                  qnorm(0.9875), 1, 0)
                ans$resid <- z$residuals/ans$scale
            }

            # VT:: 26.12.2004 
            if(p > 1)
                names(ans$coefficients)[2:p] <- xn[1:(p - 1)]
            names(ans$coefficients)[1] <- xn[p]
            
            ans$crit <- sum(z$residuals^2)
            if (intercept) {
                s1 <- sum(z$residuals^2)
                center <- mean(y)
                sh <- sum((y - center)^2)
                ans$rsquared <- 1 - (s1/sh)
            }
            else {
                s1 <- sum(z$residuals^2)
                sh <- sum(y^2)
                ans$rsquared <- 1 - (s1/sh)
            }
            if (ans$rsquared > 1) {
                ans$rsquared <- 1
            }
            if (ans$rsquared < 0) {
                ans$rsquared <- 0
            }
            ans$residuals <- rep(NA, length(na.y))
            ans$residuals[ok] <- z$residuals
            ans$lts.wt <- matrix(NA, length(na.y))
            ans$lts.wt[ok] <- weights
            ans$intercept <- intercept
            ans$method <- paste("Least Squares Regression.")
            if (abs(s0) < 1e-07) 
                ans$method <- paste(ans$method, "\nAn exact fit was found!")
            if (mcd) {
                # vt:: changed name of the function
                # mcd <- cov.mcd.default(X, print.it = FALSE, alpha = 1)
                mcd <- covMcd(X, print.it = FALSE, alpha = 1)
                if(-(determinant(mcd$cov, log = TRUE)$modulus[1])/p > 50) {
                  ans$RD[1] <- "singularity"
                }else {
                  ans$RD <- rep(NA, length(na.y))
                  ans$RD[ok] <- sqrt(mahalanobis(X, mcd$center, mcd$cov))
                  names(ans$RD) <- rownames
                }
            }
            names(ans$residuals) <- rownames
            names(ans$lts.wt) <- rownames
            names(ans$scale) <- names(ans$raw.scale) <- yn
            names(ans$rsquared) <- yn
            names(ans$crit) <- yn
            ans$X <- x
            ans$Y <- y          # VT:: 01.09.2004 - add y to the result object
            if (length(rownames)) 
                dimnames(ans$X)[[1]] <- rownames[ok]
            else {
                xx <- seq(1, length(na.x))
                dimnames(ans$X) <- list(NULL, NULL)
                dimnames(ans$X)[[1]] <- xx[ok]
            }
            ans$fitted.values <- rep(NA, length(na.y))
            ans$fitted.values[ok] <- fitted
            names(ans$fitted.values) <- rownames
            if (qr.out) 
                ans$qr <- z$qr
            class(ans) <- "lts"
            attr(ans, "call") <- sys.call()

#cat("+++++ B - alpha == 1...Ready and return\n")    
            return(ans)
        }
    }
    
    
    coefs <- rep(NA, p)
    names(coefs) <- xn
    if(qr.out)
        qrx <- qr(x)
    else 
        qrx <- qr(x)[c("rank", "pivot")]

    rk <- qrx$rank
    if (rk < p) {
        stop("x is singular")
    }
    else 
        piv <- 1:p

    if (!length(alpha)) {
        alpha <- 1/2
        quan <- quan.f(alpha, n, rk)
    }else {
        if (alpha < 1/2) 
            stop("alpha is out of range!")
        quan <- quan.f(alpha, n, rk)
    }
    y <- as.matrix(y)
    x1 <- matrix(0, ncol = p + 1, nrow = n)
    x1 <- cbind(x, y)
    x1 <- as.matrix(x1)
    storage.mode(x1) <- "double"
    datt <- matrix(0, ncol = p + 1, nrow = n)
    storage.mode(datt) <- "double"
    nvad <- p + 1
    inbest <- matrix(10000, nrow = quan, ncol = 1)
    storage.mode(inbest) <- "integer"
    objfct <- 0

    interc <- ifelse(intercept, 1, 0)
    intadjust <- ifelse(adjust, 1, 0)

    storage.mode(interc) <- "integer"
    storage.mode(seed) <- "integer"
    z <- .Fortran("rfltsreg", 
                x1 = x1, 
                as.integer(n), 
                as.integer(p), 
                as.integer(quan), 
                as.integer(nsamp), 
                inbest = inbest, 
                objfct = as.double(objfct), 
                as.integer(interc), 
                as.integer(intadjust),
                as.integer(nvad), 
                datt,
                seed,
                PACKAGE="rrcov")

    # vt:: lm.fit.qr == lm.fit(...,method=qr,...)
    #  cf <- lm.fit.qr(x[z$inbest, , drop = FALSE], y[z$inbest])$coef
    cf <- lm.fit(x[z$inbest, , drop = FALSE], y[z$inbest])$coef
    ans$best <- sort(as.vector(z$inbest))
    fitted <- x %*% cf
    resid <- y - fitted
    coefs[piv] <- cf

    # VT:: 26.12.2004 
    if(p > 1)
        ans$raw.coefficients[2:p] <- coefs[1:(p - 1)]
    ans$raw.coefficients[1] <- coefs[p]

    if(p > 1)
        names(ans$raw.coefficients)[2:p] <- names(coefs)[1:(p - 1)]
    names(ans$raw.coefficients)[1] <- names(coefs)[p]
    
    ans$alpha <- alpha
    ans$quan <- quan
    correct <- correctiefactor.s(p, intercept = intercept, n, alpha)
    s0 <- sqrt((1/quan) * sum(sort(resid^2, quan)[1:quan]))
    sh0 <- s0
    s0 <- s0 * (1/sqrt(1 - ((2 * n)/(quan * (1/qnorm((quan + n)/
        (2 * n))))) * dnorm(1/(1/(qnorm((quan + n)/(2 * n))))))) * correct
    weights <- rep(NA, n)
    if (abs(s0) < 1e-07) {
        weights <- ifelse(abs(resid) <= 1e-07, 1, 0)
        ans$scale <- ans$raw.scale <- 0
        ans$coefficients <- ans$raw.coefficients
    }
    else {
        ans$raw.scale <- s0
        ans$raw.resid <- resid/ans$raw.scale
        quantiel <- qnorm(0.9875)
        weights <- ifelse(abs(resid/s0) <= quantiel, 1, 0)
                
        # vt:: weights has to be a vector instead of a matrix - 
        #      to avoid "Error in x * wtmult : non-conformable arrays"
        # 
        weights <- as.vector(weights)
        
        z1 <- lsfit(x, y, wt = weights, intercept = FALSE)

        # VT:: 26.12.2004 
        if(p > 1)
            ans$coefficients[2:p] <- z1$coef[1:(p - 1)]
        ans$coefficients[1] <- z1$coef[p]
        
        fitted <- x %*% z1$coef
        resid <- z1$residuals
        ans$scale <- sqrt(sum(weights * resid^2)/(sum(weights) - 
            1))
        if (sum(weights) == n) {
            cdelta.rew <- 1
            correct.rew <- 1
        }
        else {
            cdelta.rew <- (1/sqrt(1 - ((2 * n)/(sum(weights) * 
                (1/qnorm((sum(weights) + n)/(2 * n))))) * dnorm(1/(1/(qnorm((sum(weights) + 
                n)/(2 * n)))))))
            correct.rew <- correctiefactor.rew.s(p, intercept = intercept, 
                n, alpha)
        }
        ans$scale <- ans$scale * cdelta.rew * correct.rew
        ans$resid <- resid/ans$scale
        quantiel <- qnorm(0.9875)
        weights <- ifelse(abs(resid/ans$scale) <= quantiel, 1, 
            0)
    }
    names(ans$coefficients) <- names(ans$raw.coefficients)
    ans$lts.wt <- matrix(NA, length(na.y))
    ans$lts.wt[ok] <- weights
    ans$crit <- z$objfct
    if (intercept) {
        initcov <- matrix(0, nrow = 1, ncol = 1)
        initmean <- matrix(0, nrow = 1, ncol = 1)
        inbest <- matrix(0, nrow = quan, ncol = 1)
        weights <- matrix(0, nrow = n, ncol = 1)
        coeff <- matrix(0, nrow = 5, ncol = 1)
        adjustcov <- matrix(0, nrow = 1, ncol = 1)
        storage.mode(initcov) <- "double"
        storage.mode(initmean) <- "double"
        storage.mode(inbest) <- "integer"
        storage.mode(weights) <- "integer"
        storage.mode(coeff) <- "double"
        storage.mode(adjustcov) <- "double"
        storage.mode(y) <- "double"
        storage.mode(seed) <- "integer"
        exactfit <- 0
        k <- 1
        sh <- .Fortran("rffastmcd", 
                    as.matrix(y), 
                    as.integer(n), 
                    as.integer(k), 
                    as.integer(quan), 
                    nsamp = 0, 
                    initcovariance = initcov, 
                    initmean = initmean, 
                    inbest = inbest, 
                    mcdestimate = 0, 
                    weights = weights, 
                    as.integer(exactfit), 
                    coeff = coeff, 
                    kount = 0, 
                    adjustcov = adjustcov,
                    seed,
                    PACKAGE="rrcov")
        y <- as.vector(y)
        sh <- as.double(sh$adjustcov)
        ans$rsquared <- 1 - (sh0/sh)^2
    }
    else {
        s1 <- sum(sort(resid^2, quan)[1:quan])
        sh <- sum(sort(y^2, quan)[1:quan])
        ans$rsquared <- 1 - (s1/sh)
    }

    # VT:: 03.11.04 - consider the case when sh=0 (i.e. rsquared=NaN or -Inf)
    if(is.nan(ans$rsquared) || is.infinite(ans$rsquared)){
        ans$rsquared <- 0
    } else if (ans$rsquared > 1) {
        ans$rsquared <- 1
    } else if (ans$rsquared < 0) {
        ans$rsquared <- 0
    }

    attributes(resid) <- attributes(fitted) <- attributes(y)
    ans$residuals <- rep(NA, length(na.y))
    ans$residuals[ok] <- resid
    ans$intercept <- intercept
    ans$method <- paste("Least Trimmed Squares Robust Regression.")
    if(abs(s0) < 1e-07) 
        ans$method <- paste(ans$method, "\nAn exact fit was found!")
    if (mcd) {

# vt:: changed name of the function
#       mcd <- cov.mcd.default(X, alpha = alpha, print.it = FALSE)
        mcd <- covMcd(X, alpha = alpha, print.it = FALSE)
        if(-(determinant(mcd$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
            ans$RD[1] <- "singularity"
        }
        else {
            ans$RD <- rep(NA, length(na.y))
            ans$RD[ok] <- sqrt(mahalanobis(X, mcd$center, mcd$cov))
            names(ans$RD) <- rownames
        }
    }
    names(ans$residuals) <- rownames
    names(ans$lts.wt) <- rownames
    names(ans$scale) <- names(ans$raw.scale) <- yn
    names(ans$rsquared) <- yn
    names(ans$crit) <- yn
    ans$X <- x
    ans$Y <- y          # VT:: 01.09.2004 - add y to the result object
    if (length(rownames)) 
        dimnames(ans$X)[[1]] <- rownames[ok]
    else {
        xx <- seq(1, length(na.x))
        dimnames(ans$X) <- list(NULL, NULL)
        dimnames(ans$X)[[1]] <- xx[ok]
    }
    ans$fitted.values <- rep(NA, length(na.y))
    ans$fitted.values[ok] <- fitted
    names(ans$fitted.values) <- rownames
    if (qr.out) 
        ans$qr <- qrx
    class(ans) <- "lts"
    attr(ans, "call") <- sys.call()
    return(ans)
}

#predict.lts <- function (object, newdata, na.action = na.pass, ...)
#{
#    if (missing(newdata)) return(fitted(object))
#    ## work hard to predict NA for rows with missing data
#    Terms <- delete.response(terms(object))
#    m <- model.frame(Terms, newdata, na.action = na.action,
#                     xlev = object$xlevels)
#    if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
#    X <- model.matrix(Terms, m, contrasts = object$contrasts)
#    drop(X %*% object$coefficients)
#} 
