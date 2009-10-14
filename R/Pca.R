setMethod("getCenter", "Pca", function(obj) obj@center)
setMethod("getLoadings", "Pca", function(obj) obj@loadings)
setMethod("getEigenvalues", "Pca", function(obj) obj@eigenvalues)
setMethod("getSdev", "Pca", function(obj) sqrt(obj@eigenvalues))
setMethod("getScores", "Pca", function(obj) obj@scores)
setMethod("getPrcomp", "Pca", function(obj) {
    ret <- list(sdev=sqrt(obj@eigenvalues),
         rotation=obj@loadings,
         center=obj@center,
         scale=FALSE,
         x=obj@scores)
    class(ret) <- "prcomp"
    ret
})

##
## Follow the standard methods: show, print, plot
##
setMethod("show", "Pca", function(object) myPcaPrint(object))
setMethod("summary", "Pca", function(object, ...){
    vars <- getEigenvalues(object)
    vars <- vars/sum(vars)
    importance <- rbind("Standard deviation" = getSdev(object),
                        "Proportion of Variance" = round(vars,5),
                        "Cumulative Proportion" = round(cumsum(vars), 5))
    colnames(importance) <- colnames(getLoadings(object))
    new("SummaryPca", pcaobj=object, importance=importance)
})
setMethod("show", "SummaryPca", function(object){

    cat("\nCall:\n")
    print(object@pcaobj@call)

    digits = max(3, getOption("digits") - 3)

    cat("Importance of components:\n")
    print(object@importance, digits = digits)
    invisible(object)
})

setMethod("print", "Pca", function(x, ...) myPcaPrint(x, ...))
setMethod("predict", "Pca", function(object, ...){
    stats:::predict.prcomp(getPrcomp(object), ...)
})
setMethod("screeplot", "Pca", function(x, ...){
    stats:::screeplot.default(getPrcomp(x), ...)
})

setMethod("biplot", "Pca", function(x, scale=1, ...){
    choices = 1L:2
    lam <- sqrt(getEigenvalues(x)[choices])
    scores <- getScores(x)
    n <- NROW(scores)
    lam <- lam * sqrt(n)

    if(scale < 0 || scale > 1)
        warning("'scale' is outside [0, 1]")

    if(scale != 0)
        lam <- lam^scale
    else
        lam <- 1

    xx <- t(t(scores[, choices]) / lam)
    yy <- t(t(getLoadings(x)[, choices]) * lam)
    .biplot(xx, yy, ...)
    invisible()
})

##  The __outlier map__ (diagnostic plot, distance-distance plot)
##  visualizes the observations by plotting their orthogonal
##  distance to the robust PCA subspace versus their robust
##  distances within the PCA subspace. This allows to classify
##  the data points  into 4 types: regular observations, good
##  leverage points, bad leverage points and orthogonal outliers.
##  The outlier plot is only possible when k < r (the number of
##  selected components is less than the rank of the matrix).
##  Otherwise a __distance plot__ will be shown (distances against
##  index).
##
##  The __screeplot__ shows the eigenvalues and is helpful to select
##  the number of principal components.
##
##  The __biplot__ is plot which aims to represent both the
##  observations and variables of a matrix of multivariate data
##  on the same plot.

## VT::17.06.2008
##setMethod("plot", "Pca", function(x, y="missing",
setMethod("plot", signature(x="Pca", y="missing"), function(x, y="missing",
                                id.n.sd=3,
                                id.n.od=3,
                                ...){

    if(all(x@od > 1.E-06))
        pca.ddplot(x, id.n.sd, id.n.od, ...)
    else
        pca.distplot(x, id.n.sd, ...)
})

myPcaPrint <- function(x, print.x=FALSE, ...) {
    if(!is.null(cl <- x@call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    cat("Standard deviations:\n"); print(sqrt(getEigenvalues(x)), ...)
    cat("\nLoadings:\n");          print(getLoadings(x), ...)

    if (print.x) {
        cat("\nRotated variables:\n"); print(getScores(x), ...)
    }
    invisible(x)
}

## Internal function to calculate the score and orthogonal distances and the
##  appropriate cutoff values for identifying outlying observations
##
##  data -
##  r    - rank
##  res  - the Pca object
##
##  - cutoff for score distances: sqrt(qchisq(0.975, k)
.distances <- function(data, r, res) {

    ## compute the score distances and the corresponding cutoff value
    n <- nrow(data)
    nk <- ncol(res@scores)
    res@sd <- sqrt(mahalanobis(res@scores, rep(0, nk), diag(res@eigenvalues, ncol=nk)))
    res@cutoff.sd <- sqrt(qchisq(0.975, res@k))

    ## Compute the orthogonal distances and the corresponding cutoff value
    ##  For each point this is the norm of the difference between the
    ##  centered data and the back-transformed scores
    res@od <- apply(data - repmat(res@center, n, 1) - res@scores %*% t(res@loadings), 1, vecnorm)
    if(is.list(dimnames(res@scores))) {
        names(res@od) <- dimnames(res@scores)[[1]]
    }

    ## The orthogonal distances make sence only of the number of PCs is less than
    ##  the rank of the data matrix - otherwise set it to 0
    res@cutoff.od <- 0
    # quan <- if(!is.null(res@quan)) res@quan else n
    quan <- getQuan(res)
    if(res@k != r) {
        ms <- unimcd(res@od^(2/3), quan=quan)
        res@cutoff.od <- sqrt(qnorm(0.975, ms$tmcd, ms$smcd)^3)
    }

    ## flag the observations with 1/0 if the distances are less or equal the
    ##  corresponding  cutoff values
    res@flag <- res@sd <= res@cutoff.sd
    if(res@cutoff.od > 0)
        res@flag <- (res@flag & res@od <= res@cutoff.od)

    return(res)
}

classSVD <- function(x){
    if(!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    else if(nrow(x) <= 1)
        stop("The sample size must be greater than 1 for svd")

    n <- nrow(x)
    p <- ncol(x)

    center <- apply(x, 2, mean)
    x <- scale(x, center=TRUE, scale=FALSE)
    svd <- svd(x/sqrt(n-1))

if(FALSE){   ## to.delete
    tolerance <- if(p < 5)          1E-12
                 else if(p <= 8)    1E-14
                 else               1E-16

    rank <- sum(svd$d > tolerance)

    tol = max(dim(x)) * svd$d[1]^2 * .Machine$double.eps
    rank <- sum(svd$d > tol)
}
    rank <- rankMM(x, sv=svd$d)
    eigenvalues <- (svd$d[1:rank])^2
    loadings <- svd$v[,1:rank]
    scores <- x %*% loadings

    list(loadings=loadings,
         scores=scores,
         eigenvalues=eigenvalues,
         rank=rank,
         x=x,
         center=center)
}

kernelEVD <- function(x){
    if(!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    else if(nrow(x) <= 1)
        stop("The sample size must be greater than 1 for svd")

    n <- nrow(x)
    p <- ncol(x)

    if(n > p) classSVD(x)
    else {
        center <- apply(x, 2, mean)
        x <- scale(x, center=TRUE, scale=FALSE)
        e <- eigen(x %*% t(x)/(n-1))

        tolerance <- n * max(e$values) * .Machine$double.eps
        rank <- sum(e$values > tolerance)

        eigenvalues <- e$values[1:rank]
        loadings <- t((x/sqrt(n-1))) %*% e$vectors[,1:rank] %*% diag(1/sqrt(eigenvalues))
        scores <- x %*% loadings
        ret <- list(loadings=loadings,
                    scores=scores,
                    eigenvalues=eigenvalues,
                    rank=rank,
                    x=x,
                    center=center)
    }
}

## Distance-distance plot (or diagnostic plot, or outlier map)
## Plots score distances against orthogonal distances
pca.ddplot <- function(obj, id.n.sd=3, id.n.od=3, title, ...) {

    if(missing(title))
    {
        title <- if(inherits(obj,"PcaClassic")) "Classical PCA" else "Robust PCA"
    }


    if(all(obj@od <= 1.E-06))
        warning("PCA diagnostic plot is not defined")
    else{
        xmax <- max(max(obj@sd), obj@cutoff.sd)
        ymax <- max(max(obj@od), obj@cutoff.od)

        plot(obj@sd, obj@od, xlab="Score distance", ylab="Orthogonal distance", xlim=c(0,xmax), ylim=c(0,ymax), type="p", ...)
        abline(v=obj@cutoff.sd)
        abline(h=obj@cutoff.od)
        label.dd(obj@sd, obj@od, id.n.sd, id.n.od)
        title(title)
    }
    invisible(obj)
}

## Distance plot, plots score distances against index
pca.distplot <- function(obj, id.n=3, title, ...) {

    if(missing(title))
    {
        title <- if(inherits(obj,"PcaClassic")) "Classical PCA" else "Robust PCA"
    }

    ymax <- max(max(obj@sd), obj@cutoff.sd)
    plot(obj@sd, xlab="Index", ylab="Score distance", ylim=c(0,ymax), type="p", ...)
    abline(h=obj@cutoff.sd)
    label(1:length(obj@sd), obj@sd, id.n)
    title(title)
    invisible(obj)
}

label <- function(x, y, id.n=3){
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    if(id.n > 0) {
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-id.n+1):n]
        if(is.character(names(y)))
            lab <- names(y[ind])
        else
            lab <- ind
        text(x[ind] + xrange/50, y[ind], lab)
    }
}

label.dd <- function(x, y, id.n.sd=3, id.n.od=3){
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    if(id.n.sd > 0 && id.n.od > 0) {
        n <- length(x)
        ind.sd <- sort(x, index.return=TRUE)$ix
        ind.sd <- ind.sd[(n - id.n.sd + 1):n]
        ind.od <- sort(y, index.return=TRUE)$ix
        ind.od <- ind.od[(n - id.n.od + 1):n]
        lab <- ind.od
        if(is.character(names(y)))
            lab <- names(y[ind.od])
        text(x[ind.od] + xrange/50, y[ind.od], lab)
        lab <- ind.sd
        if(is.character(names(x)))
            lab <- names(x[ind.sd])
        text(x[ind.sd] + xrange/50, y[ind.sd], lab)
    }
}

## VT::30.09.2009 - add a parameter 'classic' to generate a default caption
##  "Robust biplot" or "Classical biplot" for a robust/classical
##  PCA object, resp.
##  --- do not use it for now ---
##
.biplot <- function(x, y, classic,
                    var.axes = TRUE,
                    col,
                    cex = rep(par("cex"), 2),
                    xlabs = NULL, ylabs = NULL,
                    expand=1,
                    xlim = NULL, ylim = NULL,
                    arrow.len = 0.1,
                    main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ...)
{
    n <- nrow(x)
    p <- nrow(y)

##    if(is.null(main))
##        main <- if(classic) "Classical biplot" else "Robust biplot"

    if(missing(xlabs))
    {
       xlabs <- dimnames(x)[[1L]]
       if(is.null(xlabs))
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if(missing(ylabs))
    {
       ylabs <- dimnames(y)[[1L]]
       if(is.null(ylabs))
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])

    if(length(cex) == 1L)
        cex <- c(cex, cex)

    pcol <- par("col")
    if(missing(col))
    {
       col <- par("col")
       if(!is.numeric(col))
            col <- match(col, palette(), nomatch=1L)
       col <- c(col, col + 1L)
    }
    else if(length(col) == 1L)
        col <- c(col, col)

    unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])

    if(missing(xlim) && missing(ylim))
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if(!is.null(main))
        op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1L],
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    text(x, xlabs, cex = cex[1L], col = col[1L], ...)
    par(new = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
     xlab = "", ylab = "", col = col[1L], ...)

##    axis(3, col = col[2L], ...)
##    axis(4, col = col[2L], ...)
##    box(col = col[1L])
    axis(3, col = pcol, ...)
    axis(4, col = pcol, ...)
    box(col = pcol)

    text(y, labels=ylabs, cex = cex[2L], col = col[2L], ...)
    if(var.axes)
    arrows(0, 0, y[,1L] * 0.8, y[,2L] * 0.8, col = col[2L], length=arrow.len)
    invisible()
}
