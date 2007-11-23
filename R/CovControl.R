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
                           mrob=rrcov:::.mrobTau, 
                           vrob=rrcov:::.vrobGK
                           )
{
    new("CovControlOgk", niter=niter, beta=beta, mrob=mrob, vrob=vrob)        
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
