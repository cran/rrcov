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
