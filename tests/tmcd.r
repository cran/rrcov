dodata <- function(nrep=1, time=FALSE, short=FALSE, full=TRUE, method = c("FASTMCD","MASS")){
##@bdescr
## Test the function covMcd() on the literature datasets: 
##
## Call covMcd() for all regression datasets available in rrcov and print:
##  - execution time (if time == TRUE)
##  - objective fucntion
##  - best subsample found (if short == false)
##  - outliers identified (with cutoff 0.975) (if short == false)
##  - estimated center and covarinance matrix if full == TRUE)
## 
##@edescr
##
##@in  nrep              : [integer] number of repetitions to use for estimating the 
##                                   (average) execution time
##@in  time              : [boolean] whether to evaluate the execution time
##@in  short             : [boolean] whether to do short output (i.e. only the 
##                                   objective function value). If short == FALSE,
##                                   the best subsample and the identified outliers are 
##                                   printed. See also the parameter full below
##@in  full              : [boolean] whether to print the estimated cente and covariance matrix 
##@in  method            : [character] select a method: one of (FASTMCD, MASS) 

    domcd <- function(x, xname, nrep=1){ 
        if(method == "MASS"){
            mcd<-cov.mcd(x)
            quan <- as.integer((dim(x)[1] + (dim(x)[2] + 1) + 1)/2)   #default: (n+p+1)/2
        }            
        else{
            mcd<-covMcd(x, print.it=FALSE)
            quan <- as.integer(mcd$quan)
        }
        crit <- log(mcd$crit)
        if(time){
           xtime <- system.time(dorep(x, nrep, method))[1]/nrep
            xres <- sprintf("%3d %3d %3d %12.6f %10.3f\n", dim(x)[1], dim(x)[2], quan, crit, xtime)
        }
        else{
            xres <- sprintf("%3d %3d %3d %12.6f\n", dim(x)[1], dim(x)[2], quan, crit)
        }
        lpad<-lname-nchar(xname)
        cat(pad.right(xname,lpad), xres)

        if(!short){
            cat("Best subsample: \n")
            print(mcd$best)
        
            ibad <- which(mcd$mcd.wt==0)
            names(ibad) <- NULL
            nbad <- length(ibad)
            cat("Outliers: ",nbad,"\n")
            if(nbad > 0)
                print(ibad)
            if(full){
                cat("-------------\n")
                print(mcd)   
            } 
            cat("--------------------------------------------------------\n")
        }
    } 

    lname <- 20
    library(rrcov)
    method <- match.arg(method)
    if(method == "MASS")
        library(MASS)


    data(heart)
    data(phosphor)
    data(stackloss)
    data(coleman)
    data(salinity)
    data(wood)

    data(hbk)

    data(brain)
    data(milk)
    data(bushfire)

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")

    cat("Data Set               n   p  Half LOG(obj)        Time\n")
    cat("========================================================\n")
    domcd(heart.x,data(heart), nrep)
    domcd(phosphor.x,data(phosphor), nrep)
    domcd(stack.x,data(stackloss), nrep)
    domcd(coleman.x,data(coleman), nrep)
    domcd(salinity.x,data(salinity), nrep)
    domcd(wood.x,data(wood), nrep)
    domcd(hbk.x,data(hbk), nrep)

    domcd(brain,data(brain), nrep)
    domcd(milk,data(milk), nrep)
    domcd(bushfire,data(bushfire), nrep)
    cat("========================================================\n")
}

dogen <- function(nrep=1, eps=0.4, method=c("FASTMCD", "MASS")){

    domcd <- function(x, nrep=1){ 
        gc()
        xtime <- system.time(dorep(x, nrep, method))[1]/nrep
        cat(sprintf("%6d %3d %10.2f\n", dim(x)[1], dim(x)[2], xtime))
        xtime   
    } 

    set.seed(1234)

    library(rrcov)
    method <- match.arg(method)
    if(method == "MASS")
        library(MASS)

    ap <- c(2, 5, 10, 20, 30)
    an <- c(100, 500, 1000, 10000, 50000)


    tottime <- 0
    cat("     n   p       Time\n")
    cat("=====================\n")
    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- an[i]
            p <- ap[j]
            if(5*p <= n)
                tottime <- tottime + domcd(gendata(n, p, eps), nrep)
        } 
    }
    
    cat("=====================\n")
    cat("Total time: ", tottime*nrep, "\n")
}

dorep <- function(x, nrep=1, method=c("FASTMCD","MASS")){ 

    method <- match.arg(method)
    for(i in 1:nrep)
    if(method == "MASS")
        cov.mcd(x)
    else
        covMcd(x)
} 

#### gendata() ####
# Generates a location contaminated multivariate 
# normal sample of n observations in p dimensions
#    (1-eps)*Np(0,Ip) + eps*Np(m,Ip)
# where 
#    m = (b,b,...,b)
# Defaults: eps=0 and b=10
#
gendata <- function(n,p,eps=0,b=10){

    if(eps < 0 || eps >= 0.5)
        stop(message="eps must be in [0,0.5)")
    X <- mvrnorm(n,rep(0,p),diag(1,nrow=p,ncol=p))
    nbad <- as.integer(eps * n)
    if(nbad > 0){
        Xbad <- mvrnorm(nbad,rep(b,p),diag(1,nrow=p,ncol=p))
        xind <- sample(n,nbad)
        for(i in 1:nbad){
            X[xind[i],] <- Xbad[i,] 
        }
    }
    X
}

pad.right <- function(z, pads)
{
### Pads spaces to right of text
    padding <- paste(rep(" ", pads), collapse = "")
    paste(z, padding, sep = "")
}

whatis<-function(x){
    if(is.data.frame(x))
        cat("Type: data.frame\n")
    else if(is.matrix(x))
        cat("Type: matrix\n")
    else if(is.vector(x))
        cat("Type: vector\n")
    else
        cat("Type: don't know\n")
}

library(rrcov)
dodata()
