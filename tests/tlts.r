dodata <- function(nrep=1, time=FALSE, short=FALSE){

    dolts <- function(x, y, xname, nrep=1, time=FALSE, short=FALSE){ 
        lts<-ltsReg(x, y)
        crit <- lts$crit
        quan <- as.integer(lts$quan)
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
            print(lts$best)

            ibad <- which(lts$lts.wt == 0)
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

    cat("Data Set               n   p  Half      obj         Time\n")
    cat("========================================================\n")
    dolts(heart.x,heart.y, data(heart), nrep, time, short)
    dolts(phosphor.x,phosphor.y, data(phosphor), nrep, time, short)
    dolts(stack.x,stack.loss, data(stackloss), nrep, time, short)
    dolts(coleman.x,coleman.y, data(coleman), nrep, time, short)
    dolts(salinity.x,salinity.y, data(salinity), nrep, time, short)
    dolts(wood.x,wood.y, data(wood), nrep, time, short)
    dolts(hbk.x,hbk.y, data(hbk), nrep, time, short)

    cat("========================================================\n")
}

dorep <- function(x, y, nrep=1, MASS=FALSE){ 

    # set mcd=FALSE - we want to time only the LTS algorithm
    for(i in 1:nrep)
    if(MASS)
            ltsreg(x,y)
        else
            ltsReg(x, y, mcd = FALSE)
} 

dogen <- function(nrep=1, eps=0.4, MASS=FALSE){

    dolts <- function(x, y, nrep=1, MASS=FALSE){ 
        gc()
        xtime <- system.time(dorep(x, y, nrep, MASS))[1]/nrep
        n <- as.integer(dim(x)[1])
        p <- as.integer(dim(x)[2] + 1)
        cat(sprintf("%6d %3d %10.2f\n", n, p, xtime))
        xtime   
    } 

    set.seed(0)

    library(MASS)
    library(rrcov)
#    ap <- c(2, 5, 10, 20, 30)
#    an <- c(100, 500, 1000, 10000, 50000)
    ap <- c(2, 3, 5, 10)
    an <- c(100, 500, 1000, 10000, 50000)


    tottime <- 0
    cat("     n   p       Time\n")
    cat("=====================\n")
    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- an[i]
            p <- ap[j]
            if(5*p <= n){
                a <- gendata(n, p, eps)
                tottime <- tottime + dolts(a$x,a$y, nrep, MASS)
            }
        } 
    }
    
    cat("=====================\n")
    cat("Total time: ", tottime*nrep, "\n")
}

#### gendata() ####
# Generates a data set with bad leverage points (outliers in x-space) 
# n observations in p dimensions acording to the model:
#   yi = Xi1+Xi2+...+ei
# where ei - N(0,1) is the error term, Xi,j for j=1...p-1 - N(0,100) are 
# the non-trivial explanatory variables and xip is the intercept term.
# The outliers in the x-space are introduced by replacing eps. percent of
# xi1 by values distributed as N(100,100).
#
# Defaults: eps=0
#
gendata <- function(n,p,eps=0){

    if(eps < 0 || eps >= 0.5)
        stop(message="eps must be in [0,0.5)")

    p <- p-1
    x <- matrix(rnorm(n*(p),0,100), c(n,p))
    y <-rowSums(x) + 1 + rnorm(n, 0, 1)
    
    nbad <- as.integer(eps * n)
    xind <- sort(sample(n,nbad))
    xbad <- rnorm(nbad,100,100)
    for(i in 1:nbad){
        x[xind[i],1] <- xbad[i] 
    }
    list(x=x, y=y, xind=xind)
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
