###  rrcov : Scalable Robust Estimators with High Breakdown Point
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  You should have received a copy of the GNU General Public License
###  along with this program; if not, write to the Free Software
### Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for 
##  providing the initial code of this function.

## hidden in namespace:
quan.f <- function(alpha, n, rk) {
    ## Compute size of subsample, given alpha -- Same function for covMcd() and ltsReg()
    n2 <- floor((n+rk+1)/2)
    floor(2 * n2 - n + 2 * (n - n2) * alpha)
}

.ncomb <- function(k,n){
#   Computes the number of combinations of k elements out of n.

    if(n < k || k <= 0)
        return (0);
        
    comb <- 1.0
    for(j in 1:k){
        fact <- ((n - j + 1.0)/(k - j + 1.0))
        comb <- comb * fact
    }
    
    return (comb)
}

covMcd <- function(x, 
                   cor=FALSE, 
                   alpha=1/2, 
                   nsamp=500, 
                   seed=0, 
                   print.it=FALSE,
                   use.correction = TRUE,
                   control)
{

    ## no need of this parameter - used only for test and comparisons. if set to FALSE
    ##  no simulated correction factors will be used
    use.simulation = TRUE

    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    if(!missing(control)){
        defcontrol <- rrcov.control()      # default control
        if(alpha == defcontrol$alpha)       alpha <- control$alpha
        if(nsamp == defcontrol$nsamp)       nsamp <- control$nsamp
        if(seed == defcontrol$seed)         seed <- control$seed
        if(print.it == defcontrol$print.it) print.it <- control$print.it
        if(use.correction == defcontrol$use.correction) use.correction <- control$use.correction
    }

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
        stop(paste("Invalid number of trials nsamp=",nsamp,"!"))
           
    ## vt:: tolerance to be used for computing the mahalanobis distances (default = 1e-7)
    tol = 1e-10

    if(is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
        if(!is.numeric(x))
            stop("x is not a numeric dataframe or matrix.")
    }
    if((!is.vector(x) && !is.matrix(x)) || is.data.frame(x)) {
        if((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x,data.class) == "numeric")))
            stop("x is not a numeric dataframe or matrix.")
    }
    
    ## vt:: if the data is supplied as a data.frame, the following expressions results in an error
    ## as workaround convert the data.frame to a matrix
    if(is.data.frame(x))
        x <- as.matrix(x)
    
    else if(!is.matrix(x)){
        x <- array(x, c(length(x), 1), list(names(x), deparse(substitute(data))))
        x <- as.matrix(x)
    }
    dimn <- dimnames(x)
    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok,  , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    n <- dx[1]
    p <- dx[2]
    if(n < 2 * p)
        stop("Need at least 2*(number of variables) observations ")
    jmin <- floor((n + p + 1)/2)
    if(alpha < 1/2) {
        stop("The MCD must cover at least", jmin, "observations")
    }
    else if(alpha > 1)
        stop("alpha is out of range")
    quan <- quan.f(alpha, n, p)

    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will 
    ##   contain the correction factors (concistency and finite sample)
    ##   for the raw and reweighted estimates respectively. Set them initially to 1.
    ##   If use.correction is set to FALSE (default=TRUE), the finite sample correction
    ##   factor will not be used (neither for the raw estimates nor for the reweighted)
    raw.cnp2 <- rep(1,2)
    cnp2 <- rep(1,2)
    
    ## alpha=1: Compute the classical estimates
    if(alpha == 1) {
        mcd <- cov.wt(x)$cov
        loc <- as.vector(colMeans(x))
        obj <- determinant(mcd, log = TRUE)$modulus[1]
        if(( - obj/p) > 50) {
            ans <- list()
            ans$cov <- mcd
            dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
            ans$center <- loc
            if(length(dimn[[2]]))
                names(ans$center) <- dimn[[2]]
            ans$n.obs <- n
            ans$call <- match.call()
            ans$method <- paste("Minimum Covariance Determinant Estimator.")
            ans$method <- paste(ans$method, "\nThe classical covariance matrix is singular.")
            if(!print.it) {
                cat("The classical covariance matrix is singular.\n")
            }
            ans$alpha <- alpha
            ans$quan <- quan
            ans$raw.cov <- mcd
            dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
            ans$raw.center <- loc
            if(length(dimn[[2]]))
                names(ans$raw.center) <- dimn[[2]]
            ans$crit <- exp(obj)
            ans$mcd.wt <- rep(NA, length(na.x))
            ans$mcd.wt[ok] <- rep(1, sum(ok == TRUE))
        }
        else {
            mah <- mahalanobis(x, loc, mcd, tol=tol)    # VT:: 01.09.2004 - bug in alpha=1 
            ## (tol instead of tol.inv as parameter name)
            weights <- ifelse(mah < qchisq(0.975, p), 1, 0)
            ans <- cov.wt(x, wt = weights, cor)
            ans$cov <- sum(weights)/(sum(weights) - 1) * ans$cov    
        
            ## Consistency factor for reweighted MCD
            if(sum(weights) == n)
                cdelta.rew <- 1
            else {
                ## VT::19.3.2007 - replace this code used several times by a function MCDcons(p, alpha)
                ## - for the reweighted cov use 'sum(weights)/n' instead of alpha
                ##
                ### qdelta.rew <- qchisq(sum(weights)/n, p)
                ### cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(weights)/n)
                ### cnp2[1] <- cdelta.rew <- 1/cdeltainvers.rew
                cnp2[1] <- cdelta.rew <- MCDcons(p, sum(weights)/n)
            }
            ans$cov <- ans$cov * cdelta.rew
            ans$call <- match.call()
            ans$method <- paste("Minimum Covariance Determinant Estimator.")
            if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
                ans$method <- paste(ans$method, "\nThe reweighted MCD scatter matrix is singular.")
                if(!print.it) {
                  cat("The reweighted MCD scatter matrix is singular.\n")
                }
                ans$alpha <- alpha
                ans$quan <- quan
                ans$raw.cov <- mcd
                dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
                ans$raw.center <- loc
                if(length(dimn[[2]]))
                    names(ans$raw.center) <- dimn[[2]]
                ans$crit <- exp(obj)
                ans$mcd.wt <- rep(NA, length(na.x))
                ans$mcd.wt[ok] <- weights
                if(length(dimn[[1]]))
                    names(ans$mcd.wt) <- dimn[[1]]
                ans$wt <- NULL
                ans$X <- x
                if(length(dimn[[1]]))
                    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
                else {
                    xx <- seq(1, length(na.x))
                    dimnames(ans$X) <- list(NULL, NULL)
                    dimnames(ans$X)[[1]] <- xx[ok]
                }
                ans$method <- paste(ans$method, "\nThe minimum covariance determinant estimates based on", n, 
                    "observations \nare equal to the classical estimates.")
                if(print.it) {
                  cat(ans$method, "\n")
                }
                
                ans$raw.cnp2 <- raw.cnp2
                ans$cnp2 <- cnp2
                class(ans) <- "mcd"
                return(ans)
            }
            else {
                ans$alpha <- alpha
                ans$quan <- quan
                ans$raw.cov <- mcd
                dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
                ans$raw.center <- loc
                if(length(dimn[[2]]))
                    names(ans$raw.center) <- dimn[[2]]
                ans$crit <- exp(obj)
                mah <- mahalanobis(x, ans$center, ans$cov, tol=tol)
            }
            ans$mcd.wt <- rep(NA, length(na.x))
            ans$mcd.wt[ok] <- ifelse(mah < qchisq(0.975, p), 1, 0)
        }
        if(length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        ans$X <- x
        if(length(dimn[[1]]))
            dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
        else {
            xx <- seq(1, length(na.x))
            dimnames(ans$X) <- list(NULL, NULL)
            dimnames(ans$X)[[1]] <- xx[ok]
        }
        ans$method <- paste(ans$method, 
            "\nThe minimum covariance determinant estimates based on",
            n, "observations \nare equal to the classical estimates."
            )
        if(print.it) {
            cat(ans$method, "\n")
        }
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } ## end {alpha=1} --

    mcd <- .fastmcd(x, quan, nsamp, seed)    
    
    ## Compute the consistency correction factor for the raw MCD 
    ##  (see calfa in Croux and Haesbroeck)
    ## VT::19.3.2007 
    ###    qalpha <- qchisq(quan/n, p)
    ###    calphainvers <- pgamma(qalpha/2, p/2 + 1)/(quan/n)
    ###    raw.cnp2[1] <- calpha <- 1/calphainvers
    raw.cnp2[1] <- calpha <- MCDcons(p, quan/n)
    raw.cnp2[2] <- correct <- MCDcnp2(p, n, alpha, use.simulation)
    if(!use.correction)         # do not use finite sample correction factor
        raw.cnp2[2] <- correct <- 1.0
        
    if(p == 1) {
    ## The number of variables is 1 - compute univariate location and scale estimates
        scale <- sqrt(calpha) * as.double(mcd$initcovariance) * sqrt(correct)
    
        center <- as.double(mcd$initmean)
        if(abs(scale - 0) < 1e-07) {
            ans <- list()
            ## VT:: 22.12.04 - ans$cov must be a matrix. For example a subsequent 
            ##   call to determinant() will raise an error if it is a double.
            ##   The same for ans$raw.cov - see below
            
            ans$cov <- matrix(0)
            names(ans$cov) <- dimn[[2]][1]
            ans$center <- center
            names(ans$center) <- dimn[[2]][1]
            ans$n.obs <- n
            ans$call <- match.call()    
            ans$method <- paste("Univariate location and scale estimation.\nMore than",
                quan, "of the observations are identical.")
            ans$alpha <- alpha
            ans$quan <- quan
            ans$raw.cov <- matrix(0)
            names(ans$raw.cov) <- dimn[[2]][1]
            ans$raw.center <- center
            names(ans$raw.center) <- dimn[[2]][1]
            ans$crit <- 0  
            ans$mcd.wt <- rep(NA, length(na.x))
            ans$mcd.wt[ok] <- as.vector(ifelse(abs(x - center) <    1e-07, 1, 0))
            if(length(dimn[[1]]))
                names(ans$mcd.wt) <- dimn[[1]]
            if(print.it) {
                cat(ans$method, "\n")
            }
            ans$X <- x
            if(length(dimn[[1]]))
                dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
            else {
                xx <- seq(1, length(na.x))
                dimnames(ans$X) <- list(NULL, NULL)
                dimnames(ans$X)[[1]] <- xx[ok]
            }
            
            ans$raw.cnp2 <- raw.cnp2
            ans$cnp2 <- cnp2
            class(ans) <- "mcd"
            return(ans)
        }
        
        ## Compute the weights for the raw MCD in case p=1
        quantiel <- qchisq(0.975,p)
        weights <- ifelse(((x - center)/scale)^2  < quantiel, 1, 0)
        ans <- cov.wt(x, wt = weights, cor = cor)
        ans$cov <- sum(weights)/(sum(weights) - 1) * ans$cov
        
        ## Apply the correction factor for the reweighted cov
        if(sum(weights) == n){
            cdelta.rew <- 1
            correct.rew <- 1
        }
        else {
            ## VT::19.3.2007 - replace this code used several times by a function MCDcons(p, alpha)
            ## - for the reweighted cov use 'sum(weights)/n' instead of alpha
            ##
            ### qdelta.rew <- qchisq(sum(weights)/n, p)
            ### cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(weights)/n)
            ### cnp2[1] <- cdelta.rew <- 1/cdeltainvers.rew
            cnp2[1] <- cdelta.rew <- MCDcons(p, sum(weights)/n)
            cnp2[2] <- correct.rew <- MCDcnp2.rew(p, n, alpha, use.simulation)
            if(!use.correction)         # do not use finite sample correction factor
                cnp2[2] <- correct.rew <- 1.0
        }
        ans$cov <- ans$cov * cdelta.rew * correct.rew
        ans$call <- match.call()
        ans$method <- paste("Univariate location and scale estimation.")
        ans$alpha <- alpha
        ans$quan <- quan
        ans$raw.cov <- matrix(scale^2, nrow=1, ncol=1)
        names(ans$raw.cov) <- dimn[[2]][1]
        ans$raw.center <- as.vector(center)
        names(ans$raw.center) <- dimn[[2]][1]
        ans$crit <- 1/(quan - 1) * sum(sort((x - as.double(mcd$initmean))^2, partial = quan)[1:quan])
        center <- ans$center
        scale <- as.vector(sqrt(ans$cov))
        ans$mcd.wt <- rep(NA, length(na.x))
        weights <- ifelse(((x - center)/scale)^2 < qchisq(0.975, p), 1, 0)

        ans$mcd.wt[ok] <- weights
        if(length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        if(print.it) {
            cat(ans$method, "\n")
        }
        ans$X <- x
        if(length(dimn[[1]]))
            dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
        else {
            xx <- seq(1, length(na.x))
            dimnames(ans$X) <- list(NULL, NULL)
            dimnames(ans$X)[[1]] <- xx[ok]
        }
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } ## end p=1

    ## else  p >= 2 : ----------------------------------------------------------

    ## Apply correction factor to the raw estimates and use them to compute weights
    mcd$initcovariance <- calpha * mcd$initcovariance * correct
    dim(mcd$initcovariance) <- c(p, p)

    msg <- paste("Minimum Covariance Determinant Estimator")
    
    ## If not all observations are in general position, i.e. more than h observations lie on
    ## a hyperplane, the program still yields the MCD location and scatter matrix, 
    ## the latter being singular (as it should be), as well as the equation of the hyperplane.
    if(mcd$exactfit != 0) {
        dim(mcd$coeff) <- c(5, p)
        ans <- list()   
        ans$cov <- mcd$initcovariance
        dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
        ans$center <- as.vector(mcd$initmean)
        if(length(dimn[[2]]))
            names(ans$center) <- dimn[[2]]
        ans$n.obs <- n
        ans$call <- match.call()
        ans$method <- msg
        if(mcd$exactfit == -1) {
            stop("The program allows for at most ", mcd$kount, 
                " observations.")
        }
        if(mcd$exactfit == -2) {
            stop("The program allows for at most ", mcd$kount, 
                " variables.")
        }
        if(mcd$exactfit == 1) {
            ans$method <- paste(ans$method, 
                "\nThe covariance matrix of the data is singular.")
            if(!print.it) {
                cat("The covariance matrix of the data is singular.\n")
            }
        }
        if(mcd$exactfit == 2) {
            ans$method <- paste(ans$method, 
                "\nThe covariance matrix has become singular during\nthe iterations of the MCD algorithm."
                )
            if(!print.it) {
                cat("The covariance matrix has become singular during\nthe iterations of the MCD algorithm.\n"
                  )
            }
        }
        if(p == 2) {
            ans$method <- paste(ans$method, "\nThere are", mcd$
                kount, 
                "observations in the entire dataset of\n", n, 
                "observations that lie on the line with equation\n",
                round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+", 
                round(mcd$coeff[1, 2], digits = 4), 
                "(x_i2-m_2)=0 \nwith (m_1,m_2) the mean of these observations."
                )
            if(!print.it) {
                cat("There are", mcd$kount, 
                  "observations in the entire dataset of\n", n, 
                  "observations that lie on the line with equation\n",
                  round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+", 
                  round(mcd$coeff[1, 2], digits = 4), 
                  "(x_i2-m_2)=0 \nwith (m_1,m_2) the mean of these observations.\n"
                  )
            }
        }
        else if(p == 3) {
            ans$method <- paste(ans$method, "\nThere are", mcd$
                kount, 
                "observations in the entire dataset of\n", n, 
                "observations that lie on the plane with equation \n",
                round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+", 
                round(mcd$coeff[1, 2], digits = 4), "(x_i2-m_2)+", 
                round(mcd$coeff[1, 3], digits = 4), 
                "(x_i3-m_3)=0 \nwith (m_1,m_2) the mean of these observations."
                )
            if(!print.it) {
                cat("There are", mcd$kount, 
                  "observations in the entire dataset of\n", 
                  n, "observations that lie on the plane with equation \n",
                  round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+", 
                  round(mcd$coeff[1, 2], digits = 4), "(x_i2-m_2)+", 
                  round(mcd$coeff[1, 3], digits = 4), "(x_i3-m_3)=0\n",
                  "with (m_1,m_2) the mean of these observations.\n"
                  )
            }
        }
        else { ##  p > 3 -----------
            ans$method <- 
            paste(ans$method, "\nThere are", 
                mcd$kount, " observations in the entire dataset of\n", 
                n, "observations that lie on the hyperplane with equation \n",
                "a_1*(x_i1-m_1)+...+a_p*(x_ip-m_p)=0 \n",
                "with (m_1,...,m_p) the mean\n",
                "of these observations and coefficients a_i equal to: \n")
            if(!print.it) {
                cat("There are", mcd$kount, 
                  " observations in the entire dataset of\n", n,
                  "observations that lie on the hyperplane with equation \na_1*(x_i1-m_1)+...+a_p*(x_ip-m_p)=0 \nwith (m_1,...,m_p) the mean\nof these observations and coefficients a_i equal to: \n"
                  )
            }
            for(i in 1:p) {
                ans$method <- 
                paste(ans$method, round(mcd$coeff[1, i], digits = 4))
            }
            if(!print.it)
                print(round(mcd$coeff[1,  ], digits = 4))
        } # end {p > 3}
        ans$alpha <- alpha
        ans$quan <- quan
        ans$raw.cov <- mcd$initcovariance
        dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
        ans$raw.center <- as.vector(mcd$initmean)
        if(length(dimn[[2]]))
            names(ans$raw.center) <- dimn[[2]]
        ans$crit <- 0
        ans$mcd.wt <- rep(NA, length(na.x))
        ans$mcd.wt[ok] <- mcd$weights
        if(length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        ans$X <- x  
        if(length(dimn[[1]]))
            dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
        else {
            xx <- seq(1, length(na.x))
            dimnames(ans$X) <- list(xx[ok], NULL)
        }
        
        if(print.it) {
            cat(ans$method, "\n")
        }
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } #end exact fit <==>  (mcd$exactfit != 0)

    ## else ------ exactfit == 0 ----------------------------------------------
    
    mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tol)
    mcd$weights <- ifelse(mah < qchisq(0.975, p), 1, 0)

    weights <- mcd$weights
    weights <- as.vector(weights)   

    ## Compute and apply the consistency correction factor for the reweighted cov
    if(sum(weights) == n){
        cdelta.rew <- 1
        correct.rew <- 1
    }
    else {
        ## VT::19.3.2007 - replace this code used several times by a function MCDcons(p, alpha)
        ## - for the reweighted cov use 'sum(weights)/n' instead of alpha
        ##
        ### qdelta.rew <- qchisq(sum(weights)/n, p)
        ### cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(weights)/n)
        ### cnp2[1] <- cdelta.rew <- 1/cdeltainvers.rew
        cnp2[1] <- cdelta.rew <- MCDcons(p, sum(weights)/n)
        cnp2[2] <- correct.rew <- MCDcnp2.rew(p, n, alpha, use.simulation)
        if(!use.correction)         # do not use finite sample correction factor
            cnp2[2] <- correct.rew <- 1.0
    }

    ans <- cov.wt(x, wt = weights, cor)
    ans$cov <- sum(weights)/(sum(weights) - 1) * ans$cov
    ans$cov <- ans$cov * cdelta.rew * correct.rew
    ans$call <- match.call()
    ans$method <- msg
    
    ## vt:: add also the best found subsample to the result list
    ans$best <- sort(as.vector(mcd$best))

    ans$alpha <- alpha
    ans$quan <- quan
    ans$raw.cov <- mcd$initcovariance
    dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
    ans$raw.center <- as.vector(mcd$initmean)
    if(length(dimn[[2]]))
        names(ans$raw.center) <- dimn[[2]]
    ans$raw.weights <- weights
    ans$crit <- mcd$mcdestimate    
    ans$raw.mah <- mahalanobis(x,ans$raw.center,ans$raw.cov, tol = tol)

    ## Check if the reweighted scatter matrix is singular. 
    if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
        ans$method <- paste(ans$method, "\nThe reweighted MCD scatter matrix is singular.")
        if(!print.it) {
            cat("The reweighted MCD scatter matrix is singular.\n")
        }
        ans$mah <- ans$raw.mah
    }
    else {
        mah <- mahalanobis(x, ans$center, ans$cov, tol = tol)
        ans$mah <- mah
        weights<- ifelse(mah< qchisq(0.975, p), 1, 0)
    }
    ans$mcd.wt <- rep(NA, length(na.x))
    ans$mcd.wt[ok] <-weights
    if(length(dimn[[1]]))
        names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
        dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
    else {
        xx <- seq(1, length(na.x))
        dimnames(ans$X) <- list(xx[ok], NULL)
    }
    if(print.it)
        cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
}

### --- Namespace hidden (but parsed once and for all) : -------------
MCDcons <- function(p, alpha){
    ## VT::19.3.2007 - replace the code used several times by a function MCDcons(p, alpha)
    ##
    ## Compute the consistency correction factor for the MCD estimate
    ##  (see calfa in Croux and Haesbroeck)
    ##  - alpha = h/n = quan/n
    ##  - use the same function for the reweighted estimates, 
    ##      but instead of 'alpha' call with 'sum(weights)/n'
    
    qalpha <- qchisq(alpha, p)
    calphainvers <- pgamma(qalpha/2, p/2 + 1)/alpha
    1/calphainvers
}

MCDcnp2 <- function(p, n, alpha, use.simulation)
{
    if(p > 2) {
        ##                              "alfaq"            "betaq"    "qwaarden"
    coeffqpkwad875 <- matrix(c(-0.455179464070565, 1.11192541278794, 2,
                                   -0.294241208320834, 1.09649329149811, 3), ncol = 2)
    coeffqpkwad500 <- matrix(c(-1.42764571687802,  1.26263336932151, 2,
                                   -1.06141115981725,  1.28907991440387, 3), ncol = 2)

    y.500 <- log( - coeffqpkwad500[1, ] / p^coeffqpkwad500[2, ] )
    y.875 <- log( - coeffqpkwad875[1, ] / p^coeffqpkwad875[2, ] )

    A.500 <- cbind(1, - log(coeffqpkwad500[3, ] * p^2))
    A.875 <- cbind(1, - log(coeffqpkwad875[3, ] * p^2))
    coeffic.500 <- solve(A.500, y.500)
    coeffic.875 <- solve(A.875, y.875)
    fp.500.n <- 1 - exp(coeffic.500[1]) / n^coeffic.500[2]
    fp.875.n <- 1 - exp(coeffic.875[1]) / n^coeffic.875[2]
    }
    else { ## p <= 2
    if(p == 2) {
        fp.500.n <- 1 - exp( 0.673292623522027) / n^0.691365864961895
        fp.875.n <- 1 - exp( 0.446537815635445) / n^1.06690782995919
    }
    if(p == 1) {
        fp.500.n <- 1 - exp( 0.262024211897096) / n^0.604756680630497
        fp.875.n <- 1 - exp(-0.351584646688712) / n^1.01646567502486
    }
    }
    
    ## VT:18.04.2007 - use simulated correction factors for several p and n: 
    ## p in [1, 20] n in [2*p, ...]
    fp.x <- MCDcnp2sim(p, n)
    if(fp.x > 0 & use.simulation)
        fp.500.n <- 1/fp.x

    stopifnot(0.5 <= alpha, alpha <= 1)
    if(alpha <= 0.875)
    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
    else ##  0.875 < alpha <= 1
    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    return(1/fp.alpha.n)
} ## end{ MCDcnp2 }

MCDcnp2.rew <- function(p, n, alpha, use.simulation)
{
    if(p > 2) {
        ##                              "alfaq"            "betaq"    "qwaarden"
    coeffrewqpkwad875 <- matrix(c(-0.544482443573914, 1.25994483222292, 2,
                                      -0.343791072183285, 1.25159004257133, 3), ncol = 2)
    coeffrewqpkwad500 <- matrix(c(-1.02842572724793,  1.67659883081926, 2,
                                      -0.26800273450853,  1.35968562893582, 3), ncol = 2)

    y.500 <- log( - coeffrewqpkwad500[1, ] / p^ coeffrewqpkwad500[2, ] )
    y.875 <- log( - coeffrewqpkwad875[1, ] / p^ coeffrewqpkwad875[2, ] )

    A.500 <- cbind(1, - log(coeffrewqpkwad500[3, ] * p^2))
    coeffic.500 <- solve(A.500, y.500)
    A.875 <- cbind(1, - log(coeffrewqpkwad875[3, ] * p^2))
    coeffic.875 <- solve(A.875, y.875)
    fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
    fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
    }
    else {
    if(p == 2) {
        fp.500.n <- 1 - exp( 3.11101712909049 ) / n^ 1.91401056721863
        fp.875.n <- 1 - exp( 0.79473550581058 ) / n^ 1.10081930350091
    }
    if(p == 1) {
        fp.500.n <- 1 - exp( 1.11098143415027 ) / n^ 1.5182890270453
        fp.875.n <- 1 - exp( -0.66046776772861) / n^ 0.88939595831888
    }
    }

    ## VT:18.04.2007 - use simulated correction factors for several p and n: 
    ## p in [1, 20] n in [2*p, ...]
    fp.x <- MCDcnp2sim.rew(p, n)
    if(fp.x > 0 & use.simulation)
        fp.500.n <- 1/fp.x

    stopifnot(0.5 <= alpha, alpha <= 1)
    if(alpha <= 0.875)
    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
    else ##  0.875 < alpha <= 1
    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
    return(1/fp.alpha.n)
} ## end{ MCDcnp2.rew }


.fastmcd <- function(x, quan, nsamp, seed)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

#   parameters for partitioning
    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini

#   vt::03.02.2006 - added options "best" and "exact" for nsamp
    if(!missing(nsamp)){
        if(is.numeric(nsamp) && (nsamp < 0 || nsamp==0 && p > 1)){
            warning(paste("Invalid number of trials nsamp=",nsamp,"!Using default.\n"))
            nsamp <- -1
        }else if(nsamp == "exact" || nsamp == "best"){
            myk <- p
            if(n > 2*nmini-1){
                warning(paste("Options 'best' and 'exact' not allowed for n greater then ",2*nmini-1," . \nUsing nsamp=",nsamp,"\n"))
                nsamp <- -1
            }else{
                nall <- .ncomb(myk, n)
                if(nall > 5000 && nsamp == "best"){
                    nsamp <- 5000
                    warning(paste("Maximum 5000 subsets allowed for option 'best'. \nComputing 5000 subsets of size ",myk," out of ",n,"\n"))
                }else{
                    nsamp <- 0  #all subsamples
                    if(nall > 5000)
                        cat(paste("Computing all ",nall," subsets of size ",myk," out of ",n, "\n This may take a very long time!\n"))
                }
            }
        }
        
        if(!is.numeric(nsamp) || nsamp == -1){  # still not defined - set it to the default
            defcontrol <- rrcov.control()       # default control
            if(!is.numeric(nsamp))
                warning(paste("Invalid number of trials nsamp=",nsamp,"!Using default nsamp=",defcontrol$nsamp,"\n"))
            nsamp <- defcontrol$nsamp           # take the default nsamp
        }
    }

    deter <-  fit <- kount <- 0
    cutoff <- qchisq(0.975,p)
    chimed <- qchisq(0.5, p)
    
    storage.mode(x) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(quan) <- "integer"
    storage.mode(nsamp) <- "integer"
    storage.mode(seed) <- "integer"

    storage.mode(deter) <- "double"
    storage.mode(fit) <- "integer"
    storage.mode(kount) <- "integer"
    storage.mode(cutoff) <- "double"
    storage.mode(chimed) <- "double"


    initcov <- matrix(0, nrow = p * p, ncol = 1)
    adcov <- matrix(0, nrow = p * p, ncol = 1)
    initmean <- matrix(0, nrow = p, ncol = 1)
    inbest <- matrix(10000, nrow = quan, ncol = 1)
    plane <- matrix(0, nrow = 5, ncol = p)
    weights <- matrix(0, nrow = n, ncol = 1)

    storage.mode(initcov) <- "double"
    storage.mode(adcov) <- "double"
    storage.mode(initmean) <- "double"
    storage.mode(inbest) <- "integer"
    storage.mode(plane) <- "double"
    storage.mode(weights) <- "integer"
    
##   Allocate temporary storage for the fortran implementation
    temp <- matrix(0, nrow = n, ncol = 1)
    index1 <- matrix(0, nrow = n, ncol = 1)
    index2 <- matrix(0, nrow = n, ncol = 1)
    nmahad <- matrix(0, nrow = n, ncol = 1)
    ndist <- matrix(0, nrow = n, ncol = 1)
    am <- matrix(0, nrow = n, ncol = 1)
    am2 <- matrix(0, nrow = n, ncol = 1)
    slutn <- matrix(0, nrow = n, ncol = 1)

    med <- matrix(0, nrow = p, ncol = 1)
    mad <- matrix(0, nrow = p, ncol = 1)
    sd <- matrix(0, nrow = p, ncol = 1)
    means <- matrix(0, nrow = p, ncol = 1)
    bmeans <- matrix(0, nrow = p, ncol = 1)
    w <- matrix(0, nrow = p, ncol = 1)
    fv1 <- matrix(0, nrow = p, ncol = 1)
    fv2 <- matrix(0, nrow = p, ncol = 1)

    rec <- matrix(0, nrow = p+1, ncol = 1)
    sscp1 <- matrix(0, nrow = (p+1)*(p+1), ncol = 1)
    cova1 <- matrix(0, nrow = p*p, ncol = 1)
    corr1 <- matrix(0, nrow = p*p, ncol = 1)
    cinv1 <- matrix(0, nrow = p*p, ncol = 1)
    cova2 <- matrix(0, nrow = p*p, ncol = 1)
    cinv2 <- matrix(0, nrow = p*p, ncol = 1)
    z <- matrix(0, nrow = p*p, ncol = 1)

    cstock <- matrix(0, nrow = 10*p*p, ncol = 1)    #(10,nvmax2)
    mstock <- matrix(0, nrow = 10*p, ncol = 1)      #(10,nvmax)
    c1stock <- matrix(0, nrow = km10*p*p, ncol = 1) #(km10,nvmax2)
    m1stock <- matrix(0, nrow = km10*p, ncol = 1)   #(km10,nvmax)
    dath <- matrix(0, nrow = nmaxi*p, ncol = 1)     #(nmaxi,nvmax) 
    
    storage.mode(temp) <- "integer"
    storage.mode(index1) <- "integer"
    storage.mode(index2) <- "integer"
    storage.mode(nmahad) <- "double"
    storage.mode(ndist) <- "double"
    storage.mode(am) <- "double"
    storage.mode(am2) <- "double"
    storage.mode(slutn) <- "double"

    storage.mode(med) <- "double"
    storage.mode(mad) <- "double"
    storage.mode(sd) <- "double"
    storage.mode(means) <- "double"
    storage.mode(bmeans) <- "double"
    storage.mode(w) <- "double"
    storage.mode(fv1) <- "double"
    storage.mode(fv2) <- "double"

    storage.mode(rec) <- "double"
    storage.mode(sscp1) <- "double"
    storage.mode(cova1) <- "double"
    storage.mode(corr1) <- "double"
    storage.mode(cinv1) <- "double"
    storage.mode(cova2) <- "double"
    storage.mode(cinv2) <- "double"
    storage.mode(z) <- "double"

    storage.mode(cstock) <- "double"
    storage.mode(mstock) <- "double"
    storage.mode(c1stock) <- "double"
    storage.mode(m1stock) <- "double"
    storage.mode(dath) <- "double"
    
   
    mcd <- .Fortran("rffastmcd",
        x,
        n,
        p,
        quan,
        nsamp,
        initcovariance = initcov,
        initmean = initmean,
        best=inbest,
        mcdestimate = deter,
        weights = weights,
        exactfit = fit,
        coeff = plane,
        kount = kount,
        adjustcov = adcov,
        seed,
        temp,
        index1,
        index2,
        nmahad,
        ndist,
        am,
        am2,
        slutn,
        med,
        mad,
        sd,
        means,
        bmeans,
        w,
        fv1,
        fv2,
        rec,
        sscp1, 
        cova1, 
        corr1, 
        cinv1, 
        cova2, 
        cinv2, 
        z,
        cstock, 
        mstock, 
        c1stock, 
        m1stock, 
        dath,
        cutoff,
        chimed,
        PACKAGE="rrcov")  

    return (list(initcovariance=mcd$initcovariance,
            initmean=mcd$initmean,
            best=mcd$best,
            mcdestimate = mcd$mcdestimate,
            weights = mcd$weights,
            exactfit = mcd$exactfit,
            coeff = mcd$coef,
            kount = mcd$kount,
            adjustcov = mcd$adjustcov))
}

##
## VT:18.04.2007 - use simulated correction factors for several p and n and alpha=0.5. 
##  p in [1, 20] n in [2*p, ...]
##  see the modifications in MCDcnp2() and MCDcnp2.rew
##
##  VT::11.05.2007 - reduce the usage of the simulated correction factors to only those that 
##  are definitvely wrong (negative or very large). This is done by: 
##      a) reducing p.max
##      b) reducing n.max
##  NB: In general, "wrong" are the factors for the reweighted matrix, but whenever a simulated
##      value for the reweighted is used, the corresponding simulated must be used for the raw too.
##
MCDcnp2sim.p.min <- 1
MCDcnp2sim.p.max <- 9      # was 20 
MCDcnp2sim.ncol  <- 20     # This is the number of column in the matrices
##                 p=    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
MCDcnp2sim.n.min    <- c(1,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40)
MCDcnp2sim.n.max   <- c( 2,  6, 10, 13, 16, 18, 20, 20, 20, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60)
##MCDcnp2sim.n.max <- c(22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60)     ## these are the right (simulated) values for n.max
MCDcnp2sim.mat <- matrix(c(1, 3.075819, 1.515999, 2.156169, 1.480742, 1.765485, 1.460206, 1.603707, 1.427429, 1.504712, 1.334528, 1.48297,  1.355308, 1.383867, 1.319241, 1.36065,  1.307467, 1.365596, 1.255259, 1.352741, 1.239381, 3.15342, 1.799889, 2.258497, 1.688312, 1.906779, 1.548203, 1.724785, 1.500873, 1.573442, 1.417137, 1.540805, 1.395945, 1.472596, 1.394247, 1.377487, 1.337394, 1.369354, 1.333378, 1.3181, 1.313813, 1.315528, 2.12777, 2.718898, 1.993509, 2.220433, 1.820585, 1.97782, 1.672455, 1.770151, 1.587478, 1.685352, 1.539295, 1.584536, 1.499487, 1.50702, 1.41952, 1.449058, 1.393042, 1.432999, 1.369964, 1.400997, 1.333824, 2.950549, 2.145387, 2.382224, 1.927077, 2.032489, 1.8371, 1.877833, 1.710891, 1.756053, 1.620778, 1.657761, 1.558978, 1.56257, 1.508633, 1.534406, 1.46709, 1.468734, 1.432529, 1.455283, 1.386975, 1.417532, 2.229573, 2.494447, 2.016117, 2.190061, 1.877996, 1.978964, 1.767284, 1.836948, 1.677372, 1.743316, 1.616383, 1.655964, 1.55484, 1.594831, 1.502185, 1.543723, 1.467005, 1.491123, 1.44402, 1.446915, 1.401578, 2.580264, 2.109121, 2.240741, 1.944719, 2.043397, 1.821808, 1.89725, 1.748788, 1.786988, 1.659333, 1.697012, 1.610622, 1.616503, 1.538529, 1.562024, 1.499964, 1.529344, 1.474519, 1.483264, 1.441552, 1.434448, 2.165233, 2.320281, 2.007836, 2.086471, 1.884052, 1.950563, 1.76926, 1.843328, 1.708941, 1.741039, 1.627206, 1.644755, 1.580563, 1.593402, 1.527312, 1.568418, 1.501462, 1.502542, 1.464583, 1.467921, 1.431141, 2.340443, 2.048262, 2.161097, 1.926082, 1.995422, 1.81446, 1.853165, 1.738533, 1.784456, 1.679444, 1.696463, 1.612931, 1.629483, 1.548186, 1.580026, 1.52198, 1.531111, 1.482914, 1.484824, 1.442726, 1.447838, 2.093386, 2.185793, 1.948989, 2.02804, 1.867137, 1.907732, 1.771923, 1.800413, 1.691612, 1.720603, 1.642705, 1.649769, 1.589028, 1.598955, 1.539759, 1.55096, 1.503965, 1.50703, 1.471349, 1.469791, 1.436959, 2.218315, 1.997369, 2.041128, 1.887059, 1.928524, 1.79626, 1.827538, 1.716748, 1.735696, 1.658329, 1.664211, 1.599286, 1.611511, 1.553925, 1.562637, 1.516805, 1.529894, 1.476064, 1.482474, 1.453253, 1.458467, 2.0247, 2.07899, 1.921976, 1.949376, 1.824629, 1.851671, 1.744713, 1.765647, 1.683525, 1.685592, 1.625113, 1.624961, 1.571921, 1.581223, 1.535257, 1.537464, 1.497165, 1.504879, 1.468682, 1.469319, 1.448344, 2.092315, 1.941412, 1.969843, 1.844093, 1.866133, 1.766145, 1.783829, 1.703613, 1.709714, 1.646078, 1.654264, 1.594523, 1.598488, 1.545105, 1.555356, 1.514627, 1.521353, 1.483958, 1.487677, 1.449191, 1.459721, 1.958987, 1.985144, 1.87739, 1.879643, 1.786823, 1.799642, 1.720015, 1.724688, 1.663539, 1.662997, 1.609267, 1.615124, 1.56746, 1.562026, 1.520586, 1.52503, 1.493008, 1.502496, 1.471983, 1.468546, 1.435064, 1.994706, 1.880348, 1.894254, 1.805827, 1.815965, 1.744296, 1.743389, 1.665481, 1.681644, 1.624466, 1.626109, 1.584028, 1.5818, 1.54376, 1.547237, 1.504878, 1.515087, 1.479032, 1.47936, 1.450758, 1.45073, 1.892685, 1.91087, 1.825301, 1.827176, 1.745363, 1.746115, 1.693373, 1.701692, 1.648247, 1.637112, 1.594648, 1.592013, 1.554849, 1.55013, 1.522186, 1.520901, 1.492606, 1.493072, 1.460868, 1.46733, 1.440956, 1.92771, 1.835696, 1.841979, 1.775991, 1.766092, 1.703807, 1.708791, 1.654985, 1.655917, 1.602388, 1.611867, 1.570765, 1.573368, 1.53419, 1.529033, 1.506767, 1.503596, 1.481126, 1.471806, 1.444917, 1.451682, 1.850262, 1.855034, 1.778997, 1.789995, 1.718871, 1.717326, 1.667357, 1.666291, 1.619743, 1.631475, 1.582624, 1.58766, 1.546302, 1.545063, 1.512222, 1.517888, 1.489127, 1.487271, 1.466722, 1.463618, 1.444137, 1.8709, 1.794033, 1.80121, 1.736376, 1.740201, 1.673776, 1.682541, 1.638153, 1.642294, 1.604417, 1.597721, 1.559534, 1.559108, 1.533942, 1.529348, 1.499517, 1.501586, 1.473147, 1.473031, 1.457615, 1.452348, 1.805753, 1.812952, 1.746549, 1.747222, 1.696924, 1.694957, 1.652157, 1.650568, 1.607807, 1.613666, 1.577295, 1.570712, 1.543704, 1.538272, 1.515369, 1.517113, 1.487451, 1.491593, 1.464514, 1.464658, 1.439359, 1.823222, 1.758781, 1.767358, 1.70872, 1.712926, 1.666956, 1.667838, 1.62077, 1.621445, 1.592891, 1.58549, 1.55603, 1.559042, 1.521501, 1.523342, 1.499913, 1.501937, 1.473359, 1.472522, 1.452613, 1.452448), 
                         ncol = MCDcnp2sim.ncol)    
##                 p=     1  2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
MCDcnp2sim.n.min.rew <- c(1, 4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40)
MCDcnp2sim.n.max.rew <- c(2, 6, 10, 13, 16, 18, 20, 20, 20, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60)
##MCDcnp2sim.n.max.rew <- c(22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60) ## these are the right (simulated) values for n.max
MCDcnp2sim.mat.rew <- matrix(c(1, 0.984724, 0.970109, 0.978037, 0.979202, 0.982933, 1.001461, 1.026651, 0.981233, 1.011895, 1.017499, 0.964323, 1.026574, 1.006594, 0.980194, 1.009828, 0.998083, 0.966173, 1.009942, 0.99916, 1.021521, 2.216302, 1.418526, 1.635601, 1.31402, 1.33975, 1.251798, 1.210917, 1.133114, 1.150666, 1.138732, 1.096822, 1.076489, 1.058343, 1.045746, 1.036743, 1.008929, 1.049537, 1.028148, 1.027297, 1.020578, 1.00074, 1.73511, 2.06681, 1.545905, 1.659655, 1.456835, 1.47809, 1.331966, 1.334229, 1.231218, 1.220443, 1.198143, 1.193965, 1.142156, 1.146231, 1.124661, 1.112719, 1.089973, 1.070606, 1.082681, 1.061243, 1.053191, 2.388892, 1.847626, 1.96998, 1.630723, 1.701272, 1.521008, 1.553057, 1.382168, 1.414555, 1.326982, 1.321403, 1.265207, 1.264856, 1.200418, 1.21152, 1.17531, 1.168536, 1.140586, 1.14457, 1.111392, 1.112031, 1.968153, 2.168931, 1.784373, 1.894409, 1.667912, 1.693007, 1.545176, 1.582428, 1.45319, 1.480559, 1.371611, 1.358541, 1.330235, 1.30264, 1.257518, 1.244156, 1.221907, 1.22455, 1.178965, 1.177855, 1.166319, 2.275891, 1.866587, 2.014249, 1.750567, 1.829363, 1.650019, 1.689043, 1.562539, 1.561359, 1.473378, 1.488554, 1.411097, 1.416527, 1.35117, 1.361044, 1.30205, 1.299037, 1.250265, 1.260083, 1.218665, 1.236027, 1.95771, 2.074066, 1.847385, 1.905408, 1.71393, 1.768425, 1.63908, 1.67234, 1.564992, 1.562337, 1.49229, 1.499573, 1.420813, 1.424067, 1.383947, 1.378726, 1.33062, 1.330071, 1.279404, 1.295302, 1.263947, 2.164121, 1.871024, 1.979485, 1.782417, 1.84489, 1.706023, 1.734857, 1.622782, 1.634869, 1.55196, 1.554423, 1.482325, 1.509195, 1.440726, 1.436328, 1.386335, 1.396277, 1.347939, 1.346732, 1.310242, 1.309371, 1.938822, 2.050409, 1.834863, 1.882536, 1.737494, 1.761608, 1.65742, 1.687579, 1.591863, 1.60158, 1.520982, 1.535234, 1.470649, 1.486485, 1.42892, 1.435574, 1.384132, 1.382329, 1.343281, 1.346581, 1.315111, 2.063894, 1.880094, 1.907246, 1.78278, 1.806648, 1.6952, 1.720922, 1.63084, 1.635274, 1.565423, 1.56171, 1.512015, 1.4986, 1.463903, 1.456588, 1.422856, 1.407325, 1.376724, 1.373923, 1.346464, 1.34259, 1.898389, 1.950406, 1.812053, 1.849175, 1.72649, 1.737651, 1.646719, 1.655112, 1.587601, 1.597894, 1.539877, 1.53329, 1.495054, 1.490548, 1.445249, 1.446037, 1.410272, 1.412274, 1.375797, 1.369604, 1.341232, 1.992488, 1.830452, 1.857314, 1.758686, 1.763822, 1.683215, 1.679543, 1.619269, 1.608512, 1.565, 1.562282, 1.498869, 1.51325, 1.470912, 1.464654, 1.427573, 1.439301, 1.402308, 1.391006, 1.37074, 1.367573, 1.855502, 1.891242, 1.77513, 1.790618, 1.706443, 1.713098, 1.642896, 1.636577, 1.580366, 1.581752, 1.542937, 1.531668, 1.487894, 1.492039, 1.460304, 1.449762, 1.4219, 1.420953, 1.390137, 1.388677, 1.360506, 1.908277, 1.802091, 1.806128, 1.723757, 1.727249, 1.659883, 1.670056, 1.605209, 1.611481, 1.558846, 1.551762, 1.512951, 1.511515, 1.468948, 1.476073, 1.441508, 1.434997, 1.412687, 1.406782, 1.380452, 1.375924, 1.811415, 1.822311, 1.740544, 1.739355, 1.68127, 1.685342, 1.620281, 1.622572, 1.579611, 1.570103, 1.529881, 1.530097, 1.490041, 1.4947, 1.457329, 1.456344, 1.423363, 1.428653, 1.399988, 1.390069, 1.376594, 1.837723, 1.76039, 1.771031, 1.697404, 1.690915, 1.634409, 1.63713, 1.589594, 1.586521, 1.552974, 1.545571, 1.505923, 1.512794, 1.477833, 1.477821, 1.444241, 1.44452, 1.419258, 1.421297, 1.394924, 1.389393, 1.779716, 1.781271, 1.706031, 1.71224, 1.655099, 1.654284, 1.608878, 1.605955, 1.565683, 1.565938, 1.523594, 1.531235, 1.492749, 1.486786, 1.457635, 1.461416, 1.432472, 1.430164, 1.404441, 1.400021, 1.378273, 1.798932, 1.735577, 1.727031, 1.671049, 1.677601, 1.624427, 1.617626, 1.579533, 1.579987, 1.544635, 1.538715, 1.504538, 1.50726, 1.477163, 1.477084, 1.450861, 1.444496, 1.428416, 1.422813, 1.400185, 1.39552, 1.750193, 1.752145, 1.690365, 1.692051, 1.642391, 1.63858, 1.600144, 1.596401, 1.558305, 1.555932, 1.525968, 1.522984, 1.491563, 1.492554, 1.467575, 1.45786, 1.437545, 1.430893, 1.413983, 1.409386, 1.391943, 1.762922, 1.701346, 1.704996, 1.6556, 1.655548, 1.611964, 1.615219, 1.569103, 1.571079, 1.540617, 1.541602, 1.503791, 1.50195, 1.478069, 1.47678, 1.452458, 1.451732, 1.429144, 1.426547, 1.40363, 1.402647), 
                         ncol = MCDcnp2sim.ncol)    

MCDcnp2sim <- function(p, n){
    ret <- 0
    pind <- p - MCDcnp2sim.p.min + 1
    if(p >= MCDcnp2sim.p.min && p <=  MCDcnp2sim.p.max && n >= MCDcnp2sim.n.min[pind] && n <= MCDcnp2sim.n.max[pind]){
        nind <- n - MCDcnp2sim.n.min[pind] + 1
        ret <- MCDcnp2sim.mat[nind, pind]
    }
    ret
}

MCDcnp2sim.rew <- function(p, n){
    ret <- 0
    pind <- p - MCDcnp2sim.p.min + 1
    if(p >= MCDcnp2sim.p.min && p <=  MCDcnp2sim.p.max && n >= MCDcnp2sim.n.min.rew[pind] && n <= MCDcnp2sim.n.max.rew[pind]){
        nind <- n - MCDcnp2sim.n.min.rew[pind] + 1
        ret <- MCDcnp2sim.mat.rew[nind, pind]
    }
    ret
}
