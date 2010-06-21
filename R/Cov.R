setMethod("isClassic", "Cov", function(obj) TRUE)
setMethod("isSingular", "Cov", function(obj) .isSingular(getCov(obj)))

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
    cat("-> Method: ", object@method, "\n")
    if(is.list(object@singularity))
        cat(strwrap(robustbase:::singularityMsg(object@singularity, object@n.obs)), sep ="\n")

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


setMethod("isClassic", "SummaryCov", function(obj) TRUE)
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
