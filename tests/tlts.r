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

dorep <- function(x, y, nrep=1){ 

    for(i in 1:nrep)
        ltsReg(x, y)
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
