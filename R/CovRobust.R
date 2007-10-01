setMethod("isClassic", "CovRobust", function(obj) FALSE)
setMethod("isClassic", "SummaryCovRobust", function(obj) FALSE)
##
## Follow the standard methods: show, summary, plot
##
setMethod("show", "CovRobust", function(object){
    cat("\nCall:\n")
    print(object@call)
    cat("-> Method: ", object@method, "\n")
    if(is.list(object@singularity))
        cat(strwrap(robustbase:::singularityMsg(object@singularity, object@n.obs)), sep ="\n")

    digits = max(3, getOption("digits") - 3)
    cat("\nRobust Estimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nRobust Estimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)
    invisible(object)
}) 

setMethod("summary", "CovRobust", function(object, ...){

    new("SummaryCovRobust", covobj=object, evals=eigen(object@cov)$values)
    
})


setMethod("show", "SummaryCovRobust", function(object){

    cat("\nCall:\n")
    print(object@covobj@call)
    
    digits = max(3, getOption("digits") - 3)
    cat("\nRobust Estimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nRobust Estimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)
    
    cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(getEvals(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nRobust Distances: \n")
    print.default(format(as.vector(getDistance(object)), digits = digits), print.gap = 2, quote = FALSE)
})

setMethod("plot", "CovRobust", function(x, y="missing", 
                                which=c("all", "dd", "distance", "qqchi2", "tolEllipsePlot", "screeplot"),
                                classic= FALSE,
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

    ccov <- Cov(data)
    md <- rd <- NULL
    if(!isSingular(ccov))
        md <- sqrt(getDistance(ccov))
    if(!isSingular(x))
        rd <- sqrt(getDistance(x))
    
    which <- match.arg(which)
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))
    
    ## distance-distance plot: here we need both robust and mahalanobis distances
    if((which == "all" || which == "dd") && !is.null(md) && !is.null(rd)) {
        .myddplot(md, rd, cutoff=cutoff, id.n=id.n) # distance-distance plot
    }

    ## index plot of mahalanobis distances
    if((which == "all" || which == "distance")  && !is.null(rd)) {
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        .mydistplot(rd, cutoff, id.n=id.n) # index plot of robust distances
        if(classic && !is.null(md)) {
            .mydistplot(md, cutoff, classic=TRUE, id.n=id.n) # index plot of mahalanobis distances
            par(opr)
        }
    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if((which == "all" || which == "qqchi2")  && !is.null(rd)) {
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        .qqplot(rd, p, cutoff=cutoff, id.n=id.n) # qq-plot of the robust distances versus the
                                                 # quantiles of the chi-squared distribution
        if(classic && !is.null(md)) {
            .qqplot(md, p, cutoff=cutoff, classic=TRUE, id.n=id.n)
                                                 # qq-plot of the mahalanobis distances
            par(opr)
        }
    }

    if(which == "all" || which == "tolEllipsePlot") {
        if(length(dim(data)) >= 2 && dim(data)[2] == 2){
            if(!is.null(rd)){
                if(classic &&  !is.null(md))
                    .tolellipse(rcov=x, ccov = ccov, cutoff=cutoff, id.n=id.n, tol=tol)
                else
                    .tolellipse(rcov=x, cutoff=cutoff, id.n=id.n, tol=tol)
            }
        }
        else
            warning("Warning: For tolerance ellipses the dimension must be 2!")             
    }

    if(which == "all" || which == "screeplot") {
        if(classic == TRUE)
            .myscreeplot(ccov=ccov, rcov=x)
        else
            .myscreeplot(rcov=x)
    }
}) ## end { plot("CovRobust") }
