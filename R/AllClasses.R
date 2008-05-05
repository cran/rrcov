setClass("PsiFun", representation(n = "numeric",
                                  p = "numeric",
                                  r = "numeric",
                                  alpha = "numeric",
                                  c1 = "numeric"))

setClass("PsiBwt", representation(M = "numeric"), 
                   contains="PsiFun")

## Define class unions for optional slots, e.g. for definition
##  of slots which will be computed on demand, like the 
##  mahalanobis/robust distances
setClassUnion("Uvector", c("vector", "NULL"))
setClassUnion("Umatrix", c("matrix", "NULL"))
setClassUnion("Ulist", c("list", "NULL"))

## This is a virtual base class for control objects. Each robust
##  method like CovMest, CovOgk, etc. will derive a subclass with
##  the necessary control parameters, e.g. CovControlMest will
##  contain the control parameters for CovMest.
setClass("CovControl", representation(trace="logical",
                                      tolSolve="numeric",
                                      "VIRTUAL"))

setClass("Cov", representation(call = "language",
                              cov = "matrix",
                              center = "vector",
                              n.obs = "numeric",
                              mah = "Uvector",
                              method = "character",
                              singularity = "Ulist",
                              X = "Umatrix",
                              "VIRTUAL")) 

setClass("SummaryCov", representation(covobj = "Cov",
                              evals = "vector")) 

setClass("CovClassic", contains="Cov") 

setClass("CovRobust", representation(iter="numeric",
                                     crit="numeric",
                                     wt="Uvector",
                                     "VIRTUAL"),
                    contains="Cov") 

setClass("SummaryCovRobust", representation(),
                    contains="SummaryCov") 

setClass("CovMest", representation(vt="vector"),
                    contains="CovRobust") 

setClass("CovMcd", representation(alpha = "numeric",
                                  quan = "numeric",
                                  best = "Uvector",
                                  raw.cov = "matrix",
                                  raw.center = "vector",
                                  raw.mah = "Uvector",
                                  raw.wt = "Uvector",
                                  raw.cnp2 = "numeric", 
                                  cnp2 = "numeric"),
                    contains="CovRobust") 

setClass("CovOgk", representation(raw.cov = "matrix",
                                  raw.center = "vector",
                                  raw.mah = "Uvector",
                                  raw.wt = "Uvector"),
                    contains="CovRobust") 

setClass("CovMve", representation(alpha = "numeric",
                                  quan = "numeric",
                                  best = "Uvector",
                                  raw.cov = "matrix",
                                  raw.center = "vector",
                                  raw.mah = "Uvector",
                                  raw.wt = "Uvector"),
                    contains="CovRobust") 
setClass("CovSest", representation(),
                    contains="CovRobust") 

## Control parameters for CovMcd
setClass("CovControlMcd", representation(alpha="numeric",
                                          nsamp="numeric",
                                          seed="Uvector",
                                          use.correction="logical"),
                           prototype = list(alpha=0.5,
                                          nsamp=500,
                                          seed=0,
                                          trace=FALSE,
                                          tolSolve=10e-14,
                                          use.correction=TRUE),
                           contains="CovControl") 
                    
## Control parameters for CovMest
setClass("CovControlMest", representation(r="numeric",
                                          arp="numeric",
                                          eps="numeric",
                                          maxiter="numeric"),
                           prototype = list(r=0.45,
                                            arp=0.05,
                                            eps=1e-3,
                                            maxiter=120,
                                            trace=FALSE,
                                            tolSolve=10e-14),
                           contains="CovControl"
) 

## Control parameters for CovOgk
##
## Computes robust univariate mu and sigmma of the vector x
##  - sigma: tau scale Yohai and Zamar (1988) - a truncated 
##      standard deviation
##  - mu: weighted mean
##
## Returns a vector of length two with the calculated mu and sigma
##
.mrobTau <- function(x, c1 = 4.5, c2 = 3.0, ...)       #c2=2.36075 
{
    m0 <- median(x)                     # MED
    s0 <- median(abs(x - m0))           # MAD
    r <- abs(x-m0)/s0
    wt <- (1 - (r/c1)^2)^2
    wt <- ifelse(r <= c1, wt, 0)        # wt = weigths w(x,c1)
    m <- sum(x*wt)/sum(wt)              # mu = weighted mean
    
    r <- (x-m)/s0
    r <- r^2
    r[r > c2^2] <- c2^2                 # rho(x,c2)
    s2 <- s0^2 / length(x) * sum(r)     # sigma = tau scale (Yohai&Zamar 1988)
                                        # truncated standard deviation
    c(m, sqrt(s2))
}

##
## Compute a robust estimate of the covariance of two random 
##  variables x1 and x2. 
## Use the estimate defined by Gnanadesikan and Kettenring (1972):
##     cov(X,Y)=1/4 * (sigma(X+Y)^2 - sigma(X-Y)^2)
##  where sigma is a robust univariate scale.
##  As sigma is used the tau scale estimate defined above - mrobTau()
##
.vrobGK <- function(x1, x2, ...)
{
  (rrcov:::.mrobTau(x1+x2, ...)[2]^2 - rrcov:::.mrobTau(x1-x2, ...)[2]^2)/4.0
}

setClass("CovControlOgk", representation(niter="numeric",
                                         beta="numeric",
                                         mrob="function",
                                         vrob="function"),
                           prototype = list(niter=2,
                                            beta=0.90,
                                            mrob=rrcov:::.mrobTau,
                                            vrob=rrcov:::.vrobGK,
                                            trace=FALSE,
                                            tolSolve=10e-14),
                           contains="CovControl"
) 
## Control parameters for CovMve
setClass("CovControlMve", representation(alpha="numeric",
                                          nsamp="numeric",
                                          seed="Uvector"),
                           prototype = list(alpha=0.5,
                                          nsamp=500,
                                          seed=0,
                                          trace=FALSE,
                                          tolSolve=10e-14),
                           contains="CovControl") 
## Control parameters for CovSest
setClass("CovControlSest", representation(bdp="numeric",
                                          nsamp="numeric",
                                          seed="Uvector",
                                          algo="character"),
                           prototype = list(bdp=0.5,
                                            nsamp=500,
                                            seed=NULL,
                                            trace=FALSE,
                                            tolSolve=10e-14,
                                            algo="sfast"),
                           contains="CovControl") 
                    

###################### ROBPCAPcaHubert ####################################
setClass("Pca", representation(call = "language",
                              center = "vector",
                              loadings = "matrix",
                              eigenvalues = "vector",
                              scores = "matrix",
                              k = "numeric",
                              sd = "Uvector",
                              od = "Uvector",
                              cutoff.sd = "numeric",
                              cutoff.od = "numeric",
                              flag = "Uvector",
                              n.obs = "numeric",
                              "VIRTUAL")) 

setClass("PcaClassic", contains="Pca") 

setClass("PcaRobust", representation("VIRTUAL"),
                    contains="Pca") 

setClass("PcaHubert", representation(alpha = "numeric",
                                  quan = "numeric"),
                    contains="PcaRobust") 
setClass("PcaLocantore", representation(),
                    contains="PcaRobust") 
setClass("PcaCov", representation(),
                    contains="PcaRobust") 
