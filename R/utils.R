.isSingular <- function(mat)
{
##    return( - (determinant(mat, log = TRUE)$modulus[1] - 0)/ncol(mat) > 50)
    p <- ncol(mat)
    if(!is.qr(mat)) 
        mat <- qr(mat)
    return(mat$rank < p)
}
