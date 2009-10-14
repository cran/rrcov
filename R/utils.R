.isSingular <- function(mat)
{
##    return( - (determinant(mat, log = TRUE)$modulus[1] - 0)/ncol(mat) > 50)
    p <- ncol(mat)
    if(!is.qr(mat))
        mat <- qr(mat)
    return(mat$rank < p)
}

.check_vars_numeric <- function(mf)
{
    ## we need to test just the columns which are actually used.
    mt <- attr(mf, "terms")
    mterms <- attr(mt, "factors")
    mterms <- rownames(mterms)[apply(mterms, 1, any)]
    any(sapply(mterms, function(x) is.factor(mf[,x]) || !is.numeric(mf[,x])))
}

vecnorm <- function(x, p=2) sum(x^p)^(1/p)

## Several Matlab-like utility functions ======================================================
## Return the square root of a symetric positive definite matrix
sqrtm <- function(a){
    ##
    ## [E D] = eig(A); sqrtm(A) = E * sqrt(D) * E'
    ##
    if(!is.matrix(a) || ncol(a) != nrow(a))
        stop("The matrix A must be a square matrix\n")

    ee <- eigen(a)
    if(any(ee$values < 0)) {
        stop("The matrix A must be positive definite.")
    }
    ee$vectors %*% diag(sqrt(ee$values)) %*% t(ee$vectors)
}

## Return an n by p matrix of ones
ones <- function(n=1, p=1){
    matrix(1, nrow=n, ncol=p)
}

## Return an n by p matrix of zeros
zeros <- function(n=1, p=1){
    matrix(0, nrow=n, ncol=p)
}

## Purpose: rank of a matrix ``as Matlab''
## ----------------------------------------------------------------------
## Arguments:  A: a numerical matrix, maybe non-square
##           tol: numerical tolerance (compared to singular values)
##            sv: vector of non-increasing singular values of A
##                (pass as argument if already known)
## ----------------------------------------------------------------------
## Author: Martin Maechler, Date:  7 Apr 2007, 16:16
rankMM <- function(A, tol = NULL, sv = svd(A,0,0)$d) {
    d <- dim(A)
    stopifnot(length(d)==2, length(sv)==min(d), diff(sv) < 0)   # must be sorted decreasingly
    if(is.null(tol))
        tol <- max(d) * .Machine$double.eps * abs(sv[1])
    else
        stopifnot(is.numeric(tol), tol >= 0)
    sum(sv >= tol)
}


####    rankM <- function(A){
####        qr(A)$rank
####    }


##  <Matlab>
##      a=[1 2 ; 3 4];
##      repmat(a,2,3)
##
##  <R>
##      a <- matrix(1:4,2,byrow=T)
##      repmat(a,2,3)

repmat <- function(a, n, m) {

    if(is.vector(a))    # we need a column matrix, not a vector, speaking in R terms
        a <- t(a)
    kronecker(matrix(1,n,m), a)
}
