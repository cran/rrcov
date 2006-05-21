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

## This is a virtual base class for control objects. Each robust
##  method like CovMest, CovOgk, etc. will derive a subclass with
##  the necessary control parameters, e.g. CovControlMest will
##  contain the control parameters for CovMest.
setClass("CovControl", representation = "VIRTUAL") 

setClass("Cov", representation(call = "language",
                              cov = "matrix",
                              center = "vector",
                              n.obs = "numeric",
                              mah = "Uvector",
                              method = "character",
                              X = "Umatrix")) 

setClass("SummaryCov", representation(covobj = "Cov",
                              evals = "vector")) 

setClass("CovRobust", representation(iter="numeric",
                                     crit="numeric",
                                     wt="vector",
                                     "VIRTUAL"),
                    contains="Cov") 

setClass("SummaryCovRobust", representation(),
                    contains="SummaryCov") 

setClass("CovMest", representation(vt="vector"),
                    contains="CovRobust") 

## Control parameters for CovMest
setClass("CovControlMest", representation(r="numeric",
                                          arp="numeric",
                                          eps="numeric",
                                          maxiter="numeric"),
                           prototype = list(r=0.45,
                                            arp=0.05,
                                            eps=1e-3,
                                            maxiter=120),
                           contains="CovControl"
) 
