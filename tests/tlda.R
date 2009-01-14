library(rrcov)
library(MASS)
dodata <- function(method) {

    options(digits = 5)
    set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("===================================================\n")

    data(hemophilia);   show(Linda(as.factor(gr)~., data=hemophilia, method=method))
    data(anorexia);     show(Linda(Treat~., data=anorexia, method=method))
    data(Pima.tr);      show(Linda(type~., data=Pima.tr, method=method))
    data(iris);         show(Linda(Species~., data=iris, method=method))
    data(crabs);        show(Linda(sp~., data=crabs, method=method))
    
    show(Linda(as.factor(gr)~., data=hemophilia))
    show(Linda(Treat~., data=anorexia))
    show(Linda(type~., data=Pima.tr))
    show(Linda(Species~., data=iris))
    show(Linda(sp~., data=crabs))
    cat("===================================================\n")
}


## -- now do it:
dodata(method="mcdA")
dodata(method="mcdB")
dodata(method="mcdC")
#dodata(method="fsa")
