dodata <- function(nrep=1, time=FALSE, short=FALSE){

    domcd <- function(x, xname, nrep=1, time=FALSE, short=FALSE){ 
        mcd<-covMcd(x, print.it=FALSE)
        quan <- mcd$quan
        crit <- log(mcd$crit)
        if(time){
           xtime <- system.time(dorep(x, nrep))[1]/nrep
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
        }
    } 

    lname <- 20
    library(rrcov)

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

    cat("Data Set               n   p  Half LOG(obj)        Time\n")
    cat("========================================================\n")
    domcd(heart.x,data(heart), nrep, time, short)
    domcd(phosphor.x,data(phosphor), nrep, time, short)
    domcd(stack.x,data(stackloss), nrep, time, short)
    domcd(coleman.x,data(coleman), nrep, time, short)
    domcd(salinity.x,data(salinity), nrep, time, short)
    domcd(wood.x,data(wood), nrep, time, short)
    domcd(hbk.x,data(hbk), nrep, time, short)

    domcd(brain,data(brain), nrep, time, short)
    domcd(milk,data(milk), nrep, time, short)
    domcd(bushfire,data(bushfire), nrep, time, short)
    cat("========================================================\n")
}

dogen <- function(nrep=1, eps=0.4){

    domcd <- function(x, nrep=1){ 
        xtime <- system.time(dorep(x, nrep))[1]/nrep
        cat(sprintf("%6d %3d %10.3f\n", dim(x)[1], dim(x)[2], xtime))
        xtime   
    } 

    library(MASS)
    library(rrcov)
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

dorep <- function(x, nrep=1){ 

    for(i in 1:nrep)
        covMcd(x, print.it=FALSE)
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
