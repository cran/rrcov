CovControlMcd <- function (alpha=0.5,
                           nsamp=500,
                           seed=NULL,
                           trace=FALSE,
                           use.correction=TRUE)
{
    new("CovControlMcd", alpha = alpha, 
                         nsamp = nsamp, 
                         seed = seed, 
                         trace = trace, 
                         use.correction = use.correction)
}

setMethod("estimate", "CovControlMcd", function(obj, x, ...) 
    CovMcd(x, control = obj, ...)
)

CovControlMest <- function (r = 0.45, 
                            arp = 0.05, 
                            eps=1e-3, 
                            maxiter=120
                           )
{
    new("CovControlMest", r = r, arp = arp, eps=eps, maxiter=maxiter)        
}

setMethod("estimate", "CovControlMest", function(obj, x, ...) 
    CovMest(x, control = obj, ...)
)

CovControlOgk <- function (niter=2, 
                           beta=0.90, 
                           mrob=NULL, 
                           vrob=rrcov:::.vrobGK,
                           smrob="scaleTau2", 
                           svrob="gk"
                           )
{
    new("CovControlOgk", niter=niter, beta=beta, mrob=mrob, vrob=vrob, smrob=smrob, svrob=svrob)
}

setMethod("estimate", "CovControlOgk", function(obj, x, ...) 
    CovOgk(x, control = obj, ...)
)

CovControlMve <- function (alpha=0.5,
                           nsamp=500,
                           seed=NULL,
                           trace=FALSE)
{
    new("CovControlMve", alpha = alpha, 
                         nsamp = nsamp, 
                         seed = seed, 
                         trace = trace)
}

setMethod("estimate", "CovControlMve", function(obj, x, ...) 
    CovMve(x, control = obj, ...)
)

CovControlSest <- function (bdp=0.5,
                            arp=0.1,
                            eps=1e-5,
                            maxiter=120,
                            nsamp=500,
                            seed=NULL,
                            trace=FALSE,
                            tolSolve=10e-14,
                            method="sfast")
{
    new("CovControlSest", bdp=bdp, arp=arp, eps=eps, maxiter=maxiter, 
        nsamp=nsamp, seed=seed, trace=trace, tolSolve=tolSolve, method=method)        
} 

setMethod("estimate", "CovControlSest", function(obj, x, ...) 
    CovSest(x, control = obj, ...)
)
