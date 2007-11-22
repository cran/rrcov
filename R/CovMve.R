CovMve <- function(x, 
                   alpha=1/2, 
                   nsamp=500, 
                   seed=NULL, 
                   trace=FALSE,
                   control)
{

    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    
    defcontrol <- CovControlMcd()       # default control
    if(!missing(control)){
        if(alpha == defcontrol@alpha)       alpha <- control@alpha
        if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
    }
    tolSolve <- defcontrol@tolSolve

    xcall <- match.call()

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ##   Options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
        stop("Invalid number of trials nsamp = ",nsamp, "!")

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
              dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    h <- h.alpha.n(alpha, n, p) # h(alpha) , the size of the subsamples
    h <- floor(n/2)
    
    if(n <= p + 1)          # ==> floor((n+p+1)/2) > n - 1  -- not Ok
        stop(if (n <= p)    # absolute barrier!
             "n <= p -- you can't be serious!"
        else "n == p+1  is too small sample size for MCD")
    ## else
    if(n < 2 * p) { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
    }
    if(h > n)
        stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
    else if(alpha > 1) 
        stop("alpha must be <= 1")
    
    
    ##  Case: alpha==1
    ##  ...
    ##  ...
    ##  Case p == 1
    
    method <- "Minimum volume ellipsoid estimator"
    mve <- .fastmve(x, h, nsamp)
    
    rcenter <- mve$center
    rcov <- mve$cov * (1 + 15/(n - p))^2    # correction as in MASS
    
    ## What about correction factors for the raw estimates?
    ##  Are they applied in the C code?
    ##
    ## Again consider p == 1
    
    ## else, i.e. p >= 2
        ## handle exact fit, i.e. not general position situtations
        ##
        ## OK, in general position and mve$cov is not singular
        ##  do reweighting

        ## FIXME: here we assume that mve$cov is not singular
        ## ----- but it could be!
        mah <- mahalanobis(x, rcenter, rcov, tol = tolSolve)
        quantiel <- qchisq(0.975, p)
        cut <- quantiel * quantile(mah, h/n)/qchisq(h/n, p)   # as in MASS

        weights <- as.numeric(mah < cut)
        sum.w <- sum(weights)

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov

        xcov <- cov.wt(x, wt = weights)
        xcov$cov <- sum.w/(sum.w - 1) * xcov$cov

        raw.mah <- mah
        raw.weights <- weights

        ## Check if the reweighted scatter matrix is singular and
        ##  compute distances and weights based on it
        if( - (determinant(xcov$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
            ## ans$singularity <- list(kind = "reweighted.MCD")
            if(trace) 
                cat("The reweighted MCD scatter matrix is singular.\n")
            mah <- raw.mah
        }
        else {
            mah <- mahalanobis(x, xcov$center, xcov$cov, tol = tolSolve)
            weights <- as.numeric(mah < quantiel)
        }

    ##
    ans <- new("CovMve", 
                call=xcall,
                iter=nsamp,
                crit=mve$scale,
                cov=xcov$cov,
                center=xcov$center,
                mah=mah,
                wt=weights,
                n.obs=n,
                X=x,
                method=method,
                best=mve$best,
                alpha=alpha,
                quan=h,
                raw.center=rcenter,
                raw.cov=rcov,
                raw.mah=raw.mah,
                raw.wt=raw.weights)
    ans
                         
}

fastmve.m <- function(x, nsamp=500, seed=99) {
    n <- nrow(x)
    p <- ncol(x)
    n2 <- floor(n/2)
    nind <- p +1
    set.seed(seed)
    tmp <- .C('r_fast_mve', 
        as.double(x), 
        as.integer(n), 
        as.integer(p), 
        as.integer(nsamp), 
        nsing = as.integer(0), 
        ctr = as.double(rep(0,p)), 
        cov = as.double(rep(0,p*p)), 
        scale = as.double(0), 
        best=as.integer(rep(0,n)), 
        as.integer(nind), 
        as.integer(n2))

    mve.cov <- matrix(tmp$cov, p, p)
    return(list(center= tmp$ctr, 
                cov=mve.cov, 
                scale=tmp$scale,
                best=tmp$best[1:n2],
                nsamp=nsamp, 
                nsing = tmp$nsing,
                quan=n2))
}


.fastmve <- function(x, h, nsamp)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    nind <- p+1
    
    set.seed(99)
    

    if(FALSE) {## ---------------------------------------
        ##   parameters for partitioning {equal to those in Fortran !!}
        kmini <- 5
        nmini <- 300
        km10 <- 10*kmini
        nmaxi <- nmini*kmini
    
        ##   Options "best" and "exact" for nsamp
        if(!missing(nsamp)) {
            if(is.numeric(nsamp) && (nsamp < 0 || nsamp == 0 && p > 1)) {
                warning("Invalid number of trials nsamp= ", nsamp, " ! Using default.\n")
                nsamp <- -1
            } else if(nsamp == "exact" || nsamp == "best") {
                myk <- p + 1 ## was 'p'; but p+1 ("nsel = nvar+1") is correct
                if(n > 2*nmini-1) {
                    warning("Options 'best' and 'exact' not allowed for n greater than ", 2*nmini-1,".\nUsing default.\n")
                    nsamp <- -1
                } else {
                    nall <- choose(n, myk)
                    if(nall > 5000 && nsamp == "best") {
                        nsamp <- 5000
                        warning("'nsamp = \"best\"' allows maximally 5000 subsets;\n",
                            "computing these subsets of size ",
                                        myk," out of ",n,"\n")
                    } else { ## "exact" or ("best"  &  nall < 5000)
                        nsamp <- 0 ## all subsamples
                        if(nall > 5000)
                        warning("Computing all ", nall, " subsets of size ",myk,
                            " out of ",n,
                            "\n This may take a very long time!\n",
                            immediate. = TRUE)
                    }
                }
            }
        
            if(!is.numeric(nsamp) || nsamp == -1) { # still not defined - set it to the default
                defCtrl <- rrcov.control() # default control
                if(!is.numeric(nsamp))
                warning("Invalid number of trials nsamp= ",nsamp,
                    " ! Using default nsamp= ",defCtrl$nsamp,"\n")
                nsamp <- defCtrl$nsamp  # take the default nsamp
            }
        }
    }   ##  ---------------------------------------------------------------

           
    tmp <- .C('r_fast_mve', 
        x = if(is.double(x)) x else as.double(x),
        as.integer(n), 
        as.integer(p), 
        as.integer(nsamp), 
        nsing = as.integer(0), 
        ctr = as.double(rep(0,p)), 
        cov = as.double(rep(0,p*p)), 
        scale = as.double(0), 
        best=as.integer(rep(0,n)), 
        as.integer(nind), 
        as.integer(h))
        
    mve.cov <- matrix(tmp$cov, p, p)
    return(list(center= tmp$ctr, 
                cov=mve.cov, 
                scale=tmp$scale,
                best=tmp$best[1:h],
                nsamp=nsamp, 
                nsing = tmp$nsing))
}
