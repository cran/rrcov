setClass("PcaCov", representation(delta = "numeric",
                                    quan = "numeric"),
                                contains="PcaRobust") 

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
    ## this is not a `standard' model-fitting function,
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

PcaCov.default <- function(x, k=0, kmax=ncol(x), corr=FALSE, cov.control = CovControlMcd(), na.action = na.fail, trace=FALSE, ...)
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
    ##
    ## verify and set the input parameters: k and kmax
    ##
    kmax <- max(min(floor(kmax), floor(n/2), rankMM(x)),1)
    if((k <- floor(k)) < 0)   
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }
    if(k != 0)
        k <- min(k, ncol(data))
    else {
        k <- min(kmax, ncol(data))
        if(trace)
            cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="") 
    }
######################################################################

    cov <- estimate(cov.control, data)
    covmat <- list(cov=getCov(cov), center=getCenter(cov), n.obs=cov@n.obs)
    if(corr)
        covmat$cor <- getCorr(cov)

    out <- princomp(cor=corr, covmat=covmat, na.action=na.action)

    scores <- predict(out, newdata=data)
    center   <- getCenter(cov)
    sdev     <- out$sdev
    scores   <- scores[, 1:k]
    loadings <- as.matrix(out$loadings)[, 1:k]
    eigenvalues  <- (sdev^2)[1:k]

######################################################################    
    names(eigenvalues) <- NULL
    if(is.list(dimnames(data)))
        rownames(scores) <- rownames(data)  # dimnames(scores)[[1]] <- dimnames(data)[[1]]
    dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaCov")
    res <- new("PcaCov", call=cl, 
                            loadings=loadings, 
                            eigenvalues=eigenvalues, 
                            center=center, 
                            scores=scores,
                            k=p,
                            n.obs=n)
               
    ## Compute distances and flags
    res <- rrcov:::.distances(x, p, res)
    return(res)
}
