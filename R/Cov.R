Cov <- function(x, unbiased = TRUE)
{
    if(is.data.frame(x))
        x <- data.matrix(x)
    else if(!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    if(n < 2 * p)
        stop("Need at least 2*(number of variables) observations ")

    call = match.call()
    method = "Classical Estimator."
    ans <- cov.wt(x)

    nobs <- nrow(x)
    if(unbiased)
        ans$cov <- (ans$cov * (nobs))/(nobs - 1)

    ## dist <- mahalanobis(x, ans$center, ans$cov)
    
##    if(corr) {
##        dimnames(ans$cor) <- dimnames(ans$cov)
##        ans$cov <- ans$cor
##    }
        
    new("Cov", call=call, cov=ans$cov, center=ans$center, n.obs=ans$n.obs, 
        method=method, X=x)
}

##
##  Follow some isXxx() methods
##

setMethod("isClassic", "Cov", function(obj) TRUE)
setMethod("isSingular", "Cov", function(obj) .isSingular(getCov(obj)))
setMethod("isClassic", "SummaryCov", function(obj) TRUE)


##  Here we want to add some getters/setters (i.e. accessor methods).
##  The first problem is a name clash: in R the accessors are usually
##  named as the member variables(slots); 
##  we need a generic with this name; if there is a function with 
##  this name, it will not work (other signature for example) - 
##  like in the example below with cov
##      slot: cov
##      generic: cov(object), method: cov(object)
##      x <- cov(object)
##      cov(object) <- x
##  An alternative would be the way it is done in Java:
##  getCov(object), setCov(object, matrix)
##
##
setMethod("getCenter", "Cov", function(obj) obj@center)
setMethod("getCov", "Cov", function(obj) obj@cov)
setMethod("getCorr", "Cov", function(obj) cov2cor(obj@cov))
setMethod("getData", "Cov", function(obj) obj@X)
setMethod("getEvals", "Cov", function(obj) 
    eigen(getCov(obj), symmetric=TRUE, only.values=TRUE)$values)

setMethod("getDistance", "Cov", function(obj){
    if(!is(obj@mah,"NULL"))
        return(obj@mah)
        
    if(is(getData(obj), "NULL"))
        stop("Cannot compute distances: no data provided")

    dd <- mahalanobis(obj@X, obj@center, obj@cov)
    eval.parent(substitute(obj@mah <- dd))
    return(dd)
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("show", "Cov", function(object){
    cat("\nCall:\n")
    print(object@call)
    
    digits = max(3, getOption("digits") - 3)
    cat("\nEstimate of Location: \n")
    print.default(format(object@center, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nEstimate of Covariance: \n")
    print.default(format(object@cov, digits = digits), print.gap = 2, quote = FALSE)
    invisible(object)
}) 

setMethod("summary", "Cov", function(object, ...){

    new("SummaryCov", covobj=object, evals=eigen(object@cov)$values)
    
})


setMethod("getCenter", "SummaryCov", function(obj) getCenter(obj@covobj))
setMethod("getCov", "SummaryCov", function(obj) getCov(obj@covobj))
setMethod("getDistance", "SummaryCov", function(obj) getDistance(obj@covobj))
setMethod("getEvals", "SummaryCov", function(obj) obj@evals)

setMethod("show", "SummaryCov", function(object){

    cat("\nCall:\n")
    print(object@covobj@call)
    
    digits = max(3, getOption("digits") - 3)
    cat("\nEstimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nEstimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)
    
    cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(getEvals(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nMahalanobis Distances: \n")
    print.default(format(as.vector(getDistance(object)), digits = digits), print.gap = 2, quote = FALSE)
})

setMethod("plot", "Cov", function(x, y="missing", 
                                which=c("all", "distance", "qqchi2", "tolEllipsePlot", "screeplot"),
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff,
                                id.n,
                                tol = 1e-7, ...)
{
    data <- getData(x)
    ##  parameters and preconditions
    if(is.vector(data) || is.matrix(data)) {
        if(!is.numeric(data))
            stop( "x is not a numeric dataframe or matrix.")
    } else if(is.data.frame(data)) {
        if(!all(sapply(data,data.class) == "numeric"))
            stop( "x is not a numeric dataframe or matrix.")
    }

    n <- dim(data)[1]
    p <- dim(data)[2]

    if(length(getCenter(x))  == 0 ||  length(getCov(x)) == 0)
        stop( "Invalid object: attributes center and cov missing!")

    if(length(getCenter(x))  != p)
        stop( "Data set and provided center have different dimensions!")

    ## Check for singularity of the cov matrix
    if(isSingular(x))
        stop("The covariance matrix is singular!")

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    md <- sqrt(getDistance(x))
    which <- match.arg(which)

    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    ## index plot of mahalanobis distances
    if(which == "all" || which == "distance") {
        .mydistplot(md, cutoff, classic=TRUE, id.n=id.n)
    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if(which == "all" || which == "qqchi2") {
        .qqplot(md, p, cutoff=cutoff, classic=TRUE, id.n=id.n)   
    }

    if(which == "all" || which == "tolEllipsePlot") {
        if(p == 2)
            .tolellipse(ccov = x, cutoff=cutoff, id.n=id.n, tol=tol)
        else
            warning("Warning: For tolerance ellipses the dimension must be 2!")
    }

    if(which == "all" || which == "screeplot") {
        .myscreeplot(ccov=x)
    }
}) ## end { plot("Cov") }
