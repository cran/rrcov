setMethod("getQuan", "PcaCov", function(obj) obj@n.obs)

##  The S3 version
PcaCov <- function (x, ...)
    UseMethod("PcaCov")

PcaCov.formula <- function (formula, data = NULL, subset, na.action, ...)
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a 'standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- PcaCov.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaCov")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

PcaCov.default <- function(x, k=0, kmax=ncol(x), cov.control = CovControlMcd(), na.action = na.fail, scale=FALSE, signflip=TRUE, trace=FALSE, ...)
{

    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    if(n < p)
        stop("'PcaCov' can only be used with more units than variables")

    ##
    ## verify and set the input parameters: k and kmax
    ##
    kmax <- max(min(floor(kmax), rankMM(x)),1)
    if((k <- floor(k)) < 0)
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }
    if(k != 0)
        k <- min(k, p)
    else {
        k <- min(kmax, p)
        if(trace)
            cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")
    }
######################################################################

    ## VT::27.08.2010: introduce 'scale' parameter; return the scale in the value object
    ##
    myscale = vector('numeric', p) + 1
    data <- scale(data, center=FALSE, scale=scale)
    mxx <- attr(data, "scaled:scale")
    if(!is.null(mxx))
        myscale <- mxx


    ## VT::30.09.2009 - add the option for classic covariance estimates - if cov.control = NULL
    covx <- if(!is.null(cov.control)) restimate(cov.control, data) else Cov(data)
    covmat <- list(cov=getCov(covx), center=getCenter(covx), n.obs=covx@n.obs)

##    if(corr)
##        covmat$cor <- getCorr(covx)
##    out <- princomp(cor=corr, covmat=covmat, na.action=na.action)

    out <- princomp(covmat=covmat, na.action=na.action)
    center   <- getCenter(covx)
    scale    <- myscale
    sdev     <- out$sdev
    loadings <- out$loadings[, 1:k, drop=FALSE]
    eigenvalues  <- (sdev^2)[1:k]

    ## VT::27.08.2010 - signflip: flip the sign of the loadings
    if(signflip)
        loadings <- .signflip(loadings)

    scores   <- scale(data, center, scale) %*% loadings
    scores   <- scores[, 1:k, drop=FALSE]

######################################################################
    names(eigenvalues) <- NULL
    if(is.list(dimnames(data)))
    {
        rownames(scores) <- rownames(data)  # dimnames(scores)[[1]] <- dimnames(data)[[1]]
    }
    dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl[[1]] <- as.name("PcaCov")
    res <- new("PcaCov", call=cl,
                            loadings=loadings,
                            eigenvalues=eigenvalues,
                            center=center,
                            scale=myscale,
                            scores=scores,
                            k=k,
                            n.obs=n)

    ## Compute distances and flags
    res <- rrcov:::.distances(x, p, res)
    return(res)
}
