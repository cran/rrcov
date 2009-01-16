setMethod("show", "Lda", function(object){

    if(!is.null(cl <- object@call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }
    
    digits = max(3, getOption("digits") - 3)
    cat("\nPrior Probabilities of Groups:\n")
    print(object@prior)
    cat("\nGroup means:\n")
    print(object@center)
    cat("\nWithin-groups Covariance Matrix:\n")
    print(object@cov)
    cat("\nLinear Coeficients:\n")
    print(object@ldf)
    cat("\nConstants:\n")
    print(object@ldfconst)

#    svd <- x$svd
#    names(svd) <- dimnames(x$scaling)[[2]]
#    if(length(svd) > 1) {
#        cat("\nProportion of trace:\n")
#        print(round(svd^2/sum(svd^2), 4), ...)
#    }
    invisible(object)
})


setMethod("predict", "Lda", function(object, newdata){

    cv <- FALSE
    if(missing(newdata))
    {
        newdata <- object@X         # use the training sample
        cv <- TRUE                  # perform cross-validation
    }
        
    x <- as.matrix(newdata)
    
    if(ncol(x) != ncol(object@center)) 
        stop("wrong number of variables")         

    ldf <- object@ldf
    ldfconst <- object@ldfconst
    ret <- .mypredict(object@prior, levels(object@grp), ldf, ldfconst, x)
    if(cv)
        ret@cv <- .confusion(object@grp, ret@classification)
   
    ret
})


.mypredict <- function(prior, lev, ldf, ldfconst, x){

    ng <- length(prior)
    nm <- names(prior)
    if(is.null(nm))
        nm <- lev

    xx <- x %*% t(ldf)
    xx <- t(t(xx) + ldfconst)
    posterior <- xx
    for(i in 1:nrow(xx)){
        tmp <- sum(exp(xx[i,]))
        for(j in 1:ncol(xx))
            posterior[i,j] <- exp(xx[i,j])/tmp
    }
    
    cl <- factor(nm[max.col(xx)], levels = lev)
    new("PredictLda", classification=cl, posterior=posterior, x = xx)
}

.AER <- function(tab)
{
    1 - sum(tab[row(tab) == col(tab)])/sum(tab)
}

.confusion <- function(actual, predicted, prior = NULL, printit=FALSE) {
    
    names <- levels(actual)
    tab <- table(actual, predicted)
    acctab <- t(apply(tab, 1, function(x) x/sum(x)))
    dimnames(acctab) <- list(Actual = names, "Predicted (cv)" = names)
    dimnames(tab) <- list(Actual = names, "Predicted (cv)" = names)
    
    if(is.null(prior)) 
    {
        cnt <- table(actual)
        prior <- cnt/sum(cnt)
    }
    else
        names(prior) <- names
        
    AER <- 1 - sum(tab[row(tab) == col(tab)])/sum(tab)

    if(printit)
    {
        prt <- as.matrix(round(c("Apparent error rate" = AER, "Prior frequency" = prior),4))
        colnames(prt) <- ""
        print(prt)
        
        cat("\nClassification table", "\n")
        print(tab)
        cat("\nConfusion matrix", "\n")
        print(round(acctab, 3))
    }
    
    invisible(tab)
}

## Internal function to perform leaving-one-out cross validation by brute force - 
##  recalculates the estimator n times, excluding each observation in turn.
##
.CV <- function(obj){

    if(!inherits(obj, "Lda"))
        stop("The object must be an Lda object")
        
    classic <- inherits(obj, "LdaClassic")
    ret <- predict(obj)
    
    X <- obj@X
    grp <- obj@grp
    ng <- length(levels(grp))
    method <- obj@method
    
    ptm <- proc.time()

    n <- nrow(X)
    p <- ncol(X)
    
    if(!classic && n*p > 500 || method == "fsa")
        warning("This could take some time!")
        
    for(i in 1:n)
    {
        cat("i=",i,"\n")
        
        ll <- if(classic)
            {
                LdaClassic(X[-i,], grouping=grp[-i]) 
            }
            else
            {
                Linda(X[-i,], grouping=grp[-i], method=method) 
            }
        
        pp <- predict(ll, newdata=t(X[i,]))
        
        ret@classification[i] <- pp@classification[1]
        ret@posterior[i,] <- pp@posterior[1,]
    }

    ret@cv <- rrcov:::.confusion(grp, ret@classification)
    
##    cat("\nElapsed time (loo): ",(proc.time() - ptm)[1],"\n")     
    ret
}

setMethod("show", "PredictLda", function(object){

    if(!is.null(object@cv))
    {
        tab <- object@cv
        acctab <- t(apply(tab, 1, function(x) x/sum(x)))
        dimnames(acctab) <- dimnames(tab)
        AER <- 1 - sum(diag(tab))/sum(tab)
        
        prt <- as.matrix(round(c("Apparent error rate" = AER),4))
        colnames(prt) <- ""
        print(prt)
        
        cat("\nClassification table", "\n")
        print(tab)
        cat("\nConfusion matrix", "\n")
        print(round(acctab, 3))    
    }
    else
        print(object@classification)
        
##    print(object@posterior)
##    print(object@x)
    invisible(object)
})

setMethod("summary", "Lda", function(object, ...){
        new("SummaryLda", ldaobj=object)
})

setMethod("show", "SummaryLda", function(object){
    show(object@ldaobj)
    invisible(object)
})
   
