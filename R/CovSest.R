## TO DO 
##
##   - 'best' and 'exact' options for nsamp, as in CovMcd
##   - start reimplementing parts in C
##
CovSest <- function(x, 
                    nsamp=500, 
                    bdp=0.5,
                    seed=NULL, 
                    trace=FALSE,
                    tolSolve=10e-14,
                    algo=c("sfast", "surreal"),
                    control)
{
    ## Compute the constant kp (used in Tbsc)
    Tbsb <- function(c, p)
    {
        ksiint <- function(c, s, p) 
        {
            (2^s)*gamma(s + p/2)*pgamma(c^2/2, s + p/2)/gamma(p/2)
        }

        y1 = ksiint(c,1,p)*3/c-ksiint(c,2,p)*3/(c^3)+ksiint(c,3,p)/(c^5);
        y2 = c*(1-pchisq(c^2,p));
        return(y1+y2)
    }

    ## Compute the tunning constant c1 for the  S-estimator with 
    ##  Tukey's biweight function given the breakdown point (bdp)
    ##  and the dimension p    
    Tbsc <- function(bdp, p)
    {
        cnew = sqrt(qchisq(1 - bdp, p))
        eps = 1e-8
        maxit = 1e3
        iter = 1
        diff = 1e6

        ## iterate until the change is less than the tollerance or
        ##  max number of iterations reached
        while(diff > eps  && iter < maxit)
        {
            cold <- cnew
            cnew <- Tbsb(cold, p)/bdp
            diff <- abs(cold - cnew)
            iter <- iter + 1
        }
        return(cnew)
    } 

    ## Analize and validate the input parameters ...
    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.

    algo <- match.arg(algo)
    if(!missing(control)){
        defcontrol <- CovControlSest()     # default control
        if(bdp == defcontrol@bdp)           bdp <- control@bdp
        if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
        if(tolSolve == defcontrol@tolSolve) tolSolve <- control@tolSolve
        if(algo == defcontrol@algo)         algo <- control@algo
    }

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    xcall <- match.call()
    n <- nrow(x)
    p <- ncol(x)
    
    if(algo == "surreal" && nsamp == 500)   # default
        nsamp = 600*p

    ## compute the constants c1 and kp
    c1 = Tbsc(bdp, p)
    kp = (c1/6) * Tbsb(c1, p)
    
    if(trace)
        cat("\nFAST-S : bdp, p, c1, kp=", bdp, p, c1, kp, "\n")

    
    mm <- if(algo == "sfast") ..fastSloc(x, nsamp=nsamp, kp=kp, cc=c1, trace=trace)
          else                ..covSURREAL(x, nsamp=nsamp, kp=kp, c1=c1, trace=trace, tol.inv=tolSolve)

    ans <- new("CovSest", 
               call = xcall,
               iter=nsamp,
               crit=mm$crit,
               cov=mm$cov,
               center=mm$center,
               n.obs=n,
               X = as.matrix(x),
               method=mm$method)
    ans
}

##
##  A fast procedure to compute an S-estimator similar to the one proposed 
##  for regression by Salibian-Barrera, M. and Yohai, V.J. (2005),
##  "A fast algorithm for S-regression estimates". This version for 
##  multivariate location/scatter is adapted from the one implemented 
##  by Kristel Joossens, K.U. Leuven, Belgium 
##  and Ella Roelant, Ghent University, Belgium.
##
## Input:
##      x      - a data matrix of size (n,p)
##      nsamp  - number of sub-samples (default=20)
##      k      - number of refining iterations in each subsample (default=2)
##      best.r - number of "best betas" to remember from the subsamples. 
##                  These will be later iterated until convergence (default=5)
##      kp, cc  - tunning constants for the  S-estimator with Tukey's biweight 
##                function given the breakdown point (bdp) and the dimension p
##
## Output: a list with components
##      center  - robust estimate of location (vector: length)
##      cov     - robust estimate of scatter (matrix: p,p)
##      crit    - value of the objective function (number)
##
..fastSloc <- function(x, nsamp, k=2, best.r=5, kp, cc, trace=FALSE) 
{
    ## NOTES:
    ##  - in the functions rho, psi, and scaledpsi=psi/u (i.e. the weight function)
    ##      is used |x| <= c1 
    ##
    ##  - function resdis() to compute the distances is used instead of 
    ##      mahalanobis() - slightly faster

    ##  The bisquare rho function:
    ##
    ##              |   x^2/2 - x^4/2*c1^2 + x^6/6*c1^4         |x| <=  c1
    ##  rho(x) =    |
    ##              |   c1^2/6                                  |x| > c1
    ##
    rho <- function(u, cc) 
    {
        w <- abs(u) <= cc
        v <- (u^2/2 * (1 - u^2/cc^2 + u^4/(3*cc^4))) * w + (1-w) * (cc^2/6)
        v
    }

    ##  The corresponding psi function: psi = rho'
    ##
    ##              |   x - 2x^3/c1^2 + x^5/c1^4         |x| <=  c1
    ##  psi(x) =    |
    ##              |   0                                |x| > c1
    ##
    ##      using ifelse is 3 times slower
    psi <- function(u, c1)
    {
        ##ifelse(abs(u) < c1, u - 2 * u^3/c1^2 + u^5/c1^4, 0)
        pp <- u - 2 * u^3/c1^2 + u^5/c1^4
        pp*(abs(u) <= c1)
    }
    
    ## weight function = psi(u)/u
    scaledpsi <- function(u, cc) 
    {
        ##ifelse(abs(xx) < c1, xx - 2 * xx^3/c1^2 + xx^5/c1^4, 0)
        pp <- (1 - (u/cc)^2)^2
        pp <- pp * cc^2/6
        pp*(abs(u) <= cc)
    }

    ## the objective function, we solve loss.S(u, s, cc) = b for "s"
    loss.S <- function(u, s, cc) 
        mean(rho(u/s, cc))

    norm <- function(x) 
        sqrt(sum(x^2))

    ## Returns square root of the mahalanobis distances of x with respect to mu and sigma
    ## Seems to be somewhat more efficient than sqrt(mahalanobis()) - by factor 1.4!
    resdis <- function(x, mu, sigma)
    {
        central <- t(x) - mu
        sqdis <- colSums(solve(sigma, central) * central)
        dis <- sqdis^(0.5)
        dis
    }

    ##  Computes Tukey's biweight objective function (scale)
    ##  (respective to the mahalanobis distances u) using the 
    ##  rho() function and the konstants kp and c1
    scaleS <- function(u, kp, c1, initial.sc=median(abs(u))/.6745) 
    {
        ## find the scale, full iterations
        max.it <- 200
        sc <- initial.sc
        i <- 0
        eps <- 1e-20
        ## magic number alert
        err <- 1
        while((i <- i+1) < max.it && err > eps) 
        {
            sc2 <- sqrt(sc^2 * mean(rho(u/sc, c1)) / kp)
            err <- abs(sc2/sc - 1)
            sc <- sc2
        }
        sc
    }

    ##
    ##  Do "k" IRWLS refining steps from "initial.mu, initial.sigma"
    ##
    ##  If "initial.scale" is present, it's used, o/w the MAD is used
    ##
    ##  k = number of refining steps
    ##  conv = 0 means "do k steps and don't check for convergence"
    ##  conv = 1 means "stop when convergence is detected, or the
    ##                 maximum number of iterations is achieved"
    ##  kp and cc = tuning constants of the equation
    ## 
    re.s <- function(x, initial.mu, initial.sigma, initial.scale, k, conv, kp, cc) 
    {
        n <- nrow(x)
        p <- ncol(x)
        rdis <- resdis(x, initial.mu, initial.sigma)
        
        if(missing(initial.scale)) 
        {
            initial.scale <- scale <- median(abs(rdis))/.6745
        } else 
        {
            scale <- initial.scale
        }
        
        ## if conv == 1 then set the max no. of iterations to 50 magic number alert!!!
        if (conv == 1) 
            k <- 50
        
        mu <- initial.mu
        sigma <- initial.sigma
        lower.bound <- median(abs(rdis))/cc
        for(i in 1:k) 
        {
            ## do one step of the iterations to solve for the scale
            scale.super.old <- scale
        
            ## lower.bound <- median(abs(rdis))/1.56
            scale <- sqrt(scale^2 * mean(rho(rdis/scale, cc)) / kp)
        
            ## now do one step of IRWLS with the "improved scale"
            weights <- scaledpsi(rdis/scale, cc)
            W <- weights %*% matrix(rep(1,p), ncol=p)
            xw <- x * W/mean(weights)
            mu.1 <- apply(xw,2,mean)
            res <- x - matrix(rep(1,n),ncol=1) %*% mu.1
            sigma.1 <- t(res) %*% ((weights %*% matrix(rep(1,p), ncol=p)) * res)
            sigma.1 <- (det(sigma.1))^(-1/p) * sigma.1
            
            ## if(det(sigma.1) < 1e-7)     # singular
            if(.isSingular(sigma.1))
            { 
                mu.1 <- initial.mu
                sigma.1 <- initial.sigma
                scale <- initial.scale
                break
            }
            if(conv == 1) 
            {
                ## check for convergence
                if(norm(mu - mu.1) / norm(mu) < 1e-20) 
                    break
                ## magic number alert!!!
            }
            
            rdis <- resdis(x,mu.1,sigma.1)
            mu <- mu.1
            sigma <- sigma.1
        }
    
        rdis <- resdis(x,mu,sigma)

        ## get the residuals from the last beta
        return(list(mu.rw = mu.1, sigma.rw=sigma.1, scale.rw = scale))
    }

################################################################################################    
    n <- nrow(x)
    p <- ncol(x)

    best.mus <- matrix(0, best.r, p)
    best.sigmas <- matrix(0,best.r*p,p)
    best.scales <- rep(1e20, best.r)
    s.worst <- 1e20
    n.ref <- 1
    
    for(i in 1:nsamp) 
    {   
        ## Generate a p+1 subsample in general position and compute 
        ##  its mu and sigma. If sigma is singular, generate another
        ##  subsample.
        ##  FIXME ---> the singularity check 'det(sigma) < 1e-7' can
        ##      in some cases give TRUE, e.g. milk although 
        ##      .isSingular(sigma) which uses qr(mat)$rank returns 
        ##      FALSE. This can result in many retrials and thus 
        ##      speeding down!
        ##
        singular <- TRUE
        while(singular) 
        {
            indices <- sample(n, p+1)
            xs <- x[indices,]
            mu <- colMeans(xs)
            sigma <- cov(xs) 
            ##singular <- det(sigma) < 1e-7
            singular <- .isSingular(sigma)
        }
        sigma <- det(sigma)^(-1/p) * sigma    

        ## Perform k steps of IRLS on the elemental set
        if(k > 0) 
        { 
            ## do the refining
            tmp <- re.s(x=x, initial.mu=mu, initial.sigma=sigma, k=k, conv=0, kp=kp, cc=cc)
            mu.rw <- tmp$mu.rw
            sigma.rw <- tmp$sigma.rw
            scale.rw <- tmp$scale.rw
            rdis.rw <- resdis(x, mu.rw, sigma.rw)
        } else 
        { 
            ## k = 0 means "no refining"
            mu.rw <- mu
            sigma.rw <- sigma
            rdis.rw <- resdis(x, mu.rw, sigma.rw)
            scale.rw <- median(abs(rdis.rw))/.6745
        }
        
        if(i > 1) 
        { 
            ## if this isn't the first iteration....
            ## check whether new mu/sigma belong to the top best results; if so keep
            ## mu and sigma with corresponding scale.
            scale.test <- loss.S(rdis.rw, s.worst, cc)
            if(scale.test < kp) 
            {
                s.best <- scaleS(rdis.rw, kp, cc, scale.rw)
                ind <- order(best.scales)[best.r]
                best.scales[ind] <- s.best
                best.mus[ind,] <- mu.rw
                bm1 <- (ind-1)*p;
                best.sigmas[(bm1+1):(bm1+p),] <- sigma.rw
                s.worst <- max(best.scales)
            }
        } else 
        { 
            ## if this is the first iteration, then this is the best solution anyway...
            best.scales[best.r] <- scaleS(rdis.rw, kp, cc, scale.rw)
            best.mus[best.r,] <- mu.rw
            bm1 <- (best.r-1)*p;
            best.sigmas[(bm1+1):(bm1+p),] <- sigma.rw
        }
    }

    ## do the complete refining step until convergence (conv=1) starting
    ## from the best subsampling candidate (possibly refined)
    super.best.scale <- 1e20
    for(i in best.r:1)
    {
        index <- (i-1)*p
        tmp <- re.s(x=x, initial.mu=best.mus[i,],
                    initial.sigma=best.sigmas[(index+1):(index+p),],
                    initial.scale=best.scales[i], 
                    k=0, conv=1, kp=kp, cc=cc)
        if(tmp$scale.rw < super.best.scale) 
        {
            super.best.scale <- tmp$scale.rw
            super.best.mu <- tmp$mu.rw
            super.best.sigma <- tmp$sigma.rw
        }   
    }

    super.best.sigma <- super.best.scale^2*super.best.sigma
    return(list(
               center=as.vector(super.best.mu),
               cov=super.best.sigma,
               crit=super.best.scale,
               method="S estimation: S-FAST"))
}

##
## Computes S-estimates of multivariate location and scatter by the Ruppert's
##  SURREAL algorithm using Tukey's biweight function
##      data    - the data: a matrix or a data.frame
##      nsamp   - number of random (p+1)-subsamples to be drawn, defaults to 600*p
##      kp, c1  - tunning constants for the  S-estimator with Tukey's biweight 
##                function given the breakdown point (bdp) and the dimension p
##
..covSURREAL <- function(data, 
                       nsamp,
                       kp,
                       c1,
                       trace=FALSE,
                       tol.inv)
{

## Compute the tunning constant c1 for the  S-estimator with 
##  Tukey's biweight function given the breakdown point (bdp)
##  and the dimension p
vcTukey<-function(bdp, p)
{
    cnew <- sqrt(qchisq(1 - bdp, p))
    maxit <- 1000
    epsilon <- 10^(-8)
    error <- 10^6
    contit <- 1
    while(error > epsilon & contit < maxit)
    {
        cold <- cnew
        kp <- vkpTukey(cold,p)
        cnew <- sqrt(6*kp/bdp)
        error <- abs(cold-cnew)
        contit <- contit+1
    }
    if(contit == maxit) 
        warning(" Maximum iterations, the last values for c1 and kp are:")
    return (list(c1 = cnew, kp = kp))
}

## Compute the constant kp (used in vcTukey)
vkpTukey<-function(c0, p)
{  
   c2 <- c0^2
   intg1 <- (p/2)*pchisq(c2, (p+2))
   intg2 <- (p+2)*p/(2*c2)*pchisq(c2, (p+4))
   intg3 <- ((p+4)*(p+2)*p)/(6*c2^2)*pchisq(c2, (p+6))
   intg4 <- (c2/6)*(1-pchisq(c2, p))
   kp <- intg1-intg2+intg3+intg4
   return(kp)
}

##  Computes Tukey's biweight objective function (scale)
##  (respective to the mahalanobis distances dx)
scaleS <- function(dx, S0, c1, kp, tol)
{
    S <- if(S0 > 0) S0 else median(dx)/qnorm(3/4)
    test <- 2*tol
    rhonew <- NA
    rhoold <- mean(rhobiweight(dx/S,c1)) - kp
    while(test >= tol)
    {
        delta <- rhoold/mean(psibiweight(dx/S,c1)*(dx/(S^2)))
        mmin <- 1
        control <- 0
        while(mmin < 10 & control != 1)
        {
            rhonew <- mean(rhobiweight(dx/(S+delta),c1))-kp
            if(abs(rhonew) < abs(rhoold))
            {
                S <- S+delta
                control <- 1
            }else
            {
                delta <- delta/2
                mmin <- mmin+1
            }   
        }
        test <- if(mmin == 10) 0 else (abs(rhoold)-abs(rhonew))/abs(rhonew)
        rhoold <- rhonew
    }
    return(abs(S))
}

rhobiweight <- function(xx, c1)
{
    lessc1 <- (xx^2)/2-(xx^4)/(2*(c1^2))+(xx^6)/(6*c1^4)
    greatc1 <- (c1^2)/6
    lessc1 * abs(xx<c1) + greatc1 * abs(xx>=c1)
}

psibiweight <- function(xx, c1)
{
    (xx - (2*(xx^3))/(c1^2)+(xx^5)/(c1^4))*(abs(xx)<c1)
}

    data <- as.matrix(data)
    n <- nrow(data)
    m <- ncol(data)

    if(missing(nsamp))
        nsamp <- 600*m

    ## TEMPORARY for testing
    ##  ignore kp and c1
    v1 <- vcTukey(0.5, m)
    c1 <- v1$c1
    kp <- v1$kp    
  
    ma <- m+1
    mb <- 1/m
    tol <- 1e-5
    Stil <- 10e10

    laux <- 1
    for(l in 1:nsamp)
    {
        ## generate a random subsample of size p+1 and compute its 
        ##  mean 'muJ' and covariance matrix 'covJ'
        singular <- TRUE
        while(singular)
        {
            xsetJ <- data[sample(1:n, ma), ]
            muJ <- apply(xsetJ, 2, mean)
            covJ <- var(xsetJ)
            singular <- .isSingular(covJ) || .isSingular(covJ <- det(covJ) ^ -mb * covJ)
        }

        if(l > ceiling(nsamp*0.2))
        {
            if(l == ceiling(nsamp*(0.5))) 
                laux <- 2
            if(l == ceiling(nsamp*(0.8))) 
                laux <- 4 

            epsilon <- runif(1)
            epsilon <- epsilon^laux
            muJ <- epsilon*muJ + (1-epsilon)*mutil
            covJ <- epsilon*covJ + (1-epsilon)*covtil
        }
        covJ <- det(covJ) ^ -mb * covJ

        md<-sqrt(mahalanobis(data, muJ, covJ, tol.inv=tol.inv))
        if(mean(rhobiweight(md/Stil, c1)) < kp)
        {
            if(Stil < 5e10)  
                Stil <- scaleS(md, Stil, c1, kp, tol)
            else 
                Stil <- scaleS(md, 0, c1, kp, tol)
            mutil <- muJ
            covtil <- covJ
            psi <- psibiweight(md, c1*Stil)
            u <- psi/md
            ubig <- diag(u)
            ux <- ubig%*%data
            muJ <- apply(ux, 2, mean)/mean(u)
            xcenter <- scale(data, center=muJ,scale=FALSE)
            covJ <- t(ubig%*%xcenter)%*%xcenter
            covJ <- det(covJ) ^ -mb * covJ

            control <- 0
            cont <- 1
            while(cont < 3 && control != 1)
            {
                cont <- cont + 1
                md <- sqrt(mahalanobis(data, muJ, covJ, tol.inv=tol.inv))
                if(mean(rhobiweight(md/Stil, c1)) < kp)
                {
                    mutil <- muJ
                    covtil <- covJ
                    control <- 1
                    if(Stil < 5e10) 
                        Stil<-scaleS(md, Stil, c1, kp, tol)
                    else 
                        Stil<-scaleS(md, 0, c1, kp, tol)
                }else
                {
                    muJ <- (muJ + mutil)/2
                    covJ <- (covJ + covtil)/2 
                    covJ <- det(covJ) ^ -mb * covJ
                }
            }
        }
    }

    covtil <- Stil^2 * covtil
    
    return (list(center = mutil,
                 cov = covtil,
                 crit = Stil,
                 method="S estimation: SURREAL")
            )
}
 
