#### -*- mode: R; kept-new-versions: 30; kept-old-versions: 20 -*-

#### MC-adjusted Outlyingness
#### ------------------------
###
### Original code from  the web site from the Antwerpen Statistics Group :
###  http://www.agoras.ua.ac.be/Robustn.htm
### which has a "MC" section and for the software links to
### ftp://ftp.win.ua.ac.be/pub/software/agoras/newfiles/mc.tar.gz
### and that contains  mcrsoft/adjoutlyingness.R

##_NEW_ (2014): moved from Antwerpen to Leuwen,
## ===> http://wis.kuleuven.be/stat/robust/software
##      has several links to 'robustbase', and S-plus code
## http://wis.kuleuven.be/stat/robust/Programs/adjusted-boxplot/adjusted-boxplot.ssc
## (copy in ../misc/Adjusted-boxplot.ssc

## MM [ FIXME ]:
## -----------

## 1)   Use  *transposed*  B[] and A[] (now called 'E') matrices   -- DONE
## 2)   use  IQR() instead of   quantile(., .75) - quantile(., .25)

##-->  but only *after* testing original code
##     ^^^^^^^^^^^^^^^^^^^^^^^^

## VT::02.08.2020 - mc() will not converge in some cases, particularly in
##  the new version of PcaHubert() in rrcov, with option for adjustment
##  for skewed data. Add a parameter mcfun which defaults to mc,
##  but can be replaced by another function, e.g. the mcVT(), given below
##  which is secured to non-convergence by tryCatch().

mcVT <- function(x)
{
    result = tryCatch({
        suppressWarnings(mc(x))
    }, error = function(e) {
        mc(x, doScale=FALSE)
    })
    result
}

adjOutlyingness <- function(x, ndir = 250, p.samp = p, clower=4, cupper=3,
                            alpha.cutoff = 0.75, coef = 1.5, qr.tol = 1e-12, keep.tol = 1e-12,
                           only.outlyingness = FALSE, maxit.mult = max(100, p), mcfun=mc, trace.lev = 0)
## Skewness-Adjusted Outlyingness
{
    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n >= 1, p >= 1, p.samp >= p, is.numeric(x))
    if (p <= n) {
        B <- matrix(0, p, ndir)
        E <- matrix(1, p, 1)
        x. <- unname(x) # for speed in subsequent subsetting and solve
        maxit <- as.integer(maxit.mult * ndir)
        ## ^^ original code had 'Inf', i.e. no iter.count check;
        ## often, maxit == ndir would suffice
        if(trace.lev >= 2) p10 <- 10 ^ max(0, min(6 - trace.lev, floor(log10(maxit))))
	i <- 1L
	it <- 0L
	while (i <= ndir  &&  (it <- it+1L) < maxit) {
            P <- x.[sample.int(n, p.samp), , drop=FALSE]
            if ((qrP <- qr(P, tol = qr.tol))$rank == p) {
                B[,i] <- solve(qrP, E, tol = qr.tol)
                if(trace.lev) cat(" it=",it,"; found direction # ", i,"\n", sep="")
		i <- i+1L
            } else if(trace.lev >= 2) {
                if(it %% p10 == 0)
                    cat(" it=",it,": rank(qr(P ..)) = ", qrP$rank, " < p = ",p,"\n", sep="")
            }
        }
	if(it == maxit) {
	    rnk.x <- qr(x, tol = qr.tol)$rank
	    if(rnk.x < p)
		stop("Matrix 'x' is not of full rank: rankM(x) = ",rnk.x," < p = ",p,
                     "\n Use fullRank(x) instead")
	    ## else
	    stop("** direction sampling iterations were not sufficient. Maybe try increasing 'maxit.mult'")
	}
        Bnorm <- sqrt(colSums(B^2))
        Nx <- mean(abs(x.)) ## so the comparison is scale-equivariant:
        keep <- Bnorm*Nx > keep.tol
        Bnormr <- Bnorm[  keep ]
        B      <-     B[, keep , drop=FALSE]

        A <- B / rep(Bnormr, each = nrow(B))
    }
    else {
        stop('More dimensions than observations: not yet implemented')
        ## MM: In LIBRA(matlab) they have it implemented:
        ##    seed=0;
        ##    nrich1=n*(n-1)/2;
        ##    ndirect=min(250,nrich1);
        ##    true = (ndirect == nrich1);
        ##    B=extradir(x,ndir,seed,true); % n*ri
        ##      ======== % Calculates ndirect directions through
        ##               % two random choosen data points from data
        ##    for i=1:size(B,1)
        ##        Bnorm(i)=norm(B(i,:),2);
        ##    end
        ##    Bnormr=Bnorm(Bnorm > 1.e-12); %ndirect*1
        ##    B=B(Bnorm > 1.e-12,:);       %ndirect*n
        ##    A=diag(1./Bnormr)*B;         %ndirect*n

    }
    Y <- x %*% A # (n x p) %*% (p, nd') == (n x nd');
    ##  nd' = ndir.final := ndir - {those not in 'keep'}

    ## Compute and sweep out the median
    med <- colMedians(Y)
    Y <- Y - rep(med, each=n)
    ## central :<==> non-adjusted  <==> "classical" outlyingness
    central <- clower == 0 && cupper == 0
    if(!central)
        ## MM: mc() could be made faster if we could tell it that med(..) = 0

        ## VT::02.08.2020 - mc() will not converge in some cases,
        ##  e.g. generate a matrix a <- matrix(rnorm(50*30), nrow=50, ncol=30)
        ##  replace the function mc() by, say mcfun() which will do tryCatch()
        ##  and if it did not converge will retry with doScale=FALSE, this
        ##  works always.

        ## tmc <- apply(Y, 2, mc) ## original Antwerpen *wrongly*: tmc <- mc(Y)
        tmc <- apply(Y, 2, mcfun) ## now we can replace the call to the medcouple function mc()
    ##                     ==

    Q13 <- apply(Y, 2, quantile, c(.25, .75), names=FALSE)
    Q1 <- Q13[1L,]; Q3 <- Q13[2L,]
    IQR <- Q3 - Q1
    ## NOTA BENE(MM): simplified definition of tup/tlo here and below
    ## 2014-10-18: "flipped" sign (which Pieter Setaert (c/o Mia H) proposed, Jul.30,2014:
    tup <- Q3 + coef*
	(if(central) IQR else IQR*exp( cupper*tmc*(tmc >= 0) + clower*tmc*(tmc < 0)))
    tlo <- Q1 - coef*
	(if(central) IQR else IQR*exp(-clower*tmc*(tmc >= 0) - cupper*tmc*(tmc < 0)))
    ## Note: all(tlo < med & med < tup) # where med = 0

    ## Instead of the loop:
    ##  for (i in 1:ndir) {
    ##      tup[i] <-  max(Y[Y[,i] < tup[i], i])
    ##      tlo[i] <- -min(Y[Y[,i] > tlo[i], i])
    ##      ## MM:  FIXED typo-bug :       ^^^ this was missing!
    ##      ## But after the fix, the function stops "working" for longley..
    ##      ## because tlo[] becomes 0 too often, YZ[.,.] = c / 0 = Inf !
    ##  }
    Yup <- Ylo <- Y
    Yup[!(Y < rep(tup, each=n))] <- -Inf
    Ylo[!(Y > rep(tlo, each=n))] <-  Inf
    tup <-  apply(Yup, 2, max) # =  max{ Y[i,] ; Y[i,] < tup[i] }
    tlo <- -apply(Ylo, 2, min) # = -min{ Y[i,] ; Y[i,] > tlo[i] }

    tY <- t(Y)
    ## Note: column-wise medians are all 0 : "x_i > m" <==> y > 0
    ## Note: this loop is pretty fast
    for (j in 1:n) { # when y = (X-med) = 0  ==> adjout = 0 rather than
	## 0 / 0 --> NaN; e.g, in  set.seed(3); adjOutlyingness(longley)
	non0 <- 0 != (y <- tY[,j]); y <- y[non0]; I <- (y > 0)
	tY[non0, j] <- abs(y) / (I*tup[non0] + (1 - I)*tlo[non0])
    }
    ## We get +Inf above for "small n"; e.g. set.seed(11); adjOutlyingness(longley)
    adjout <- apply(tY, 2, function(x) max(x[is.finite(x)]))

    if(only.outlyingness)
	adjout
    else {
	Qadj <- quantile(adjout, probs = c(1 - alpha.cutoff, alpha.cutoff))
	mcadjout <- if(cupper != 0) mc(adjout) else 0
	##			    ===
	cutoff <- Qadj[2] + coef * (Qadj[2] - Qadj[1]) *
	    (if(mcadjout > 0) exp(cupper*mcadjout) else 1)

	list(adjout = adjout, iter = it, ndir.final = sum(keep),
	     MCadjout = mcadjout, Qalph.adjout = Qadj, cutoff = cutoff,
	     nonOut = (adjout <= cutoff))
    }
}
