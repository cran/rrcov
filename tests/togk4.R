library(rrcov)
dodata <- function(nrep=1, time=FALSE, short=FALSE, full=TRUE, method = c("FASTMCD","MASS")){
    domcd <- function(x, xname, nrep=1){
        n <- dim(x)[1]
        p <- dim(x)[2]

        mcd<-CovOgk(x)
        
        xres <- sprintf("%3d %3d\n", dim(x)[1], dim(x)[2])

        lpad<-lname-nchar(xname)
        cat(pad.right(xname,lpad), xres)

        dist <- getDistance(mcd)
        quantiel <- qchisq(0.975, p)
        ibad <- which(dist >= quantiel)
        names(ibad) <- NULL
        nbad <- length(ibad)
        cat("Outliers: ",nbad,"\n")
        if(nbad > 0)
            print(ibad)
        cat("-------------\n")
        show(mcd)   
        cat("--------------------------------------------------------\n")
    } 

    lname <- 20
    library(rrcov)
    method <- match.arg(method)

    data(heart)
    data(stars)
    data(phosphor)
    data(stackloss)
    data(coleman)
    data(salinity)
    data(wood)

    data(hbk)

    data(brain)
    data(milk)
    data(bushfire)

#    data(x1000)
#    data(x5000)
    
    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")

    cat("Data Set               n   p  Half LOG(obj)        Time\n")
    cat("========================================================\n")
    domcd(heart.x,data(heart), nrep)
    domcd(stars,data(stars), nrep)
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
#    domcd(x1000$X,data(x1000), nrep)
#    domcd(x5000$X,data(x5000), nrep)
}

pad.right <- function(z, pads)
{
### Pads spaces to right of text
    padding <- paste(rep(" ", pads), collapse = "")
    paste(z, padding, sep = "")
}

dodata()
