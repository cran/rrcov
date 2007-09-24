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
