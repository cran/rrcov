test_rsquared <- function(){
    x1 <- c(2, 3, 4, 8, 12, 22, 28, 29, 33, 34, 38, 40, 41, 47, 48, 50, 51, 54, 56, 59) 
    y1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5) 
    print(ltsReg(x1,y1,alpha=0.8))
    print(ltsReg(y1,x1,alpha=0.8))
    print(ltsReg(y1,x1,alpha=0.8, intercept=FALSE))
}

dodata <- function(nrep=1, time=FALSE, short=FALSE, full=TRUE, method = c("FASTLTS","MASS")){
##@bdescr
## Test function ltsReg() on the literature datasets: 
##
## Call ltsReg() for all regression datasets available in rrcov and print:
##  - execution time (if time == TRUE)
##  - objective fucntion
##  - best subsample found (if short == false)
##  - outliers identified (with cutoff 0.975) (if short == false)
##  - estimated coeficients and scale (if full == TRUE)
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
##@in  full              : [boolean] whether to print the estimated coeficients and scale 
##@in  method            : [character] select a method: one of (FASTLTS, MASS) 

    dolts <- function(x, y, xname, nrep=1){ 
        if(method == "MASS"){
            lts <- ltsreg(x,y)
            quan <- as.integer((dim(x)[1] + (dim(x)[2] + 1) + 1)/2)   #default: (n+p+1)/2
        } else {
            lts <- ltsReg(x, y, mcd = FALSE)
            quan <- as.integer(lts$quan)
        }

        crit <- lts$crit
        if(time){
            xtime <- system.time(dorep(x, y, nrep, method))[1]/nrep
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
            if(full){
                cat("-------------\n")
                print(lts)   
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
    data(stars)
    data(phosphor)
    data(stackloss)
    data(coleman)
    data(salinity)
    data(aircraft)
    data(delivery)
    data(wood)

    data(hbk)

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    
    cat("========================================================\n")
    cat("Data Set               n   p  Half      obj         Time\n")
    cat("========================================================\n")
    dolts(heart.x,heart.y, data(heart), nrep)
    dolts(stars.x,stars.y, data(stars), nrep)
    dolts(phosphor.x,phosphor.y, data(phosphor), nrep)
    dolts(stack.x,stack.loss, data(stackloss), nrep)
    dolts(coleman.x,coleman.y, data(coleman), nrep)
    dolts(salinity.x,salinity.y, data(salinity))
    dolts(aircraft.x,aircraft.y, data(aircraft))
    dolts(delivery.x,delivery.y, data(delivery))
    dolts(wood.x,wood.y, data(wood), nrep)
    dolts(hbk.x,hbk.y, data(hbk), nrep)

    cat("========================================================\n")
}

dorep <- function(x, y, nrep=1, method=c("FASTLTS","MASS")){ 

    # set mcd=FALSE - we want to time only the LTS algorithm
    for(i in 1:nrep)
    if(method == "MASS")
#        ltsreg(x,y,control=list(psamp = NA, nsamp = "best", adjust = FALSE))
        ltsreg(x,y)
    else
        ltsReg(x, y, mcd = FALSE)
} 

dogen <- function(nrep=1, eps=0.4, method=c("FASTLTS","MASS")){

    dolts <- function(x, y, nrep=1){ 
        gc()
        xtime <- system.time(dorep(x, y, nrep, method))[1]/nrep
        n <- as.integer(dim(x)[1])
        p <- as.integer(dim(x)[2] + 1)
        cat(sprintf("%6d %3d %10.2f\n", n, p, xtime))
        xtime   
    } 

    library(rrcov)
    method <- match.arg(method)
    if(method == "MASS")
        library(MASS)

    ap <- c(2, 3, 5, 10)
    an <- c(100, 500, 1000, 10000, 20000, 30000, 50000)

    set.seed(0)

    tottime <- 0

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    
    cat("     n   p       Time\n")
    cat("=====================\n")
    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- an[i]
            p <- ap[j]
            if(5*p <= n){
                a <- gendata(n, p, eps)
                tottime <- tottime + dolts(a$x,a$y, nrep)
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
test_rsquared()
