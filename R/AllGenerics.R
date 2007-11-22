setGeneric("plot")
setGeneric("summary")

setGeneric("psi", function(obj, x) standardGeneric("psi")) 
setGeneric("wt", function(obj, x) standardGeneric("wt")) 
setGeneric("vt", function(obj, x) standardGeneric("vt")) 
setGeneric("erho", function(obj) standardGeneric("erho")) 
setGeneric("erhoLim", function(obj) standardGeneric("erhoLim")) 
setGeneric("erhoLimD", function(obj) standardGeneric("erhoLimD")) 
setGeneric("arpLim", function(obj) standardGeneric("arpLim")) 
setGeneric("csolve", function(obj) standardGeneric("csolve")) 
setGeneric("iterM", function(obj, x, t1, s, eps, maxiter) standardGeneric("iterM")) 

setGeneric("isClassic", function(obj) standardGeneric("isClassic")) 
setGeneric("isSingular", function(obj) standardGeneric("isSingular")) 

if(!isGeneric("getCenter"))
    setGeneric("getCenter", function(obj) standardGeneric("getCenter")) 
setGeneric("getCov", function(obj) standardGeneric("getCov")) 
setGeneric("getCorr", function(obj) standardGeneric("getCorr")) 
setGeneric("getData", function(obj) standardGeneric("getData")) 
setGeneric("getDistance", function(obj) standardGeneric("getDistance")) 
setGeneric("getEvals", function(obj) standardGeneric("getEvals")) 

setGeneric("estimate", function(obj, x, ...) standardGeneric("estimate")) 

if(!isGeneric("predict"))
    setGeneric("predict", function(object, ...) standardGeneric("predict")) 
    
if(!isGeneric("screeplot"))
    setGeneric("screeplot", function(x, ...) standardGeneric("screeplot")) 
    
if(!isGeneric("biplot"))
    setGeneric("biplot", function(x, ...) standardGeneric("biplot")) 
    
if(!isGeneric("getQuan"))
    setGeneric("getQuan", function(obj) standardGeneric("getQuan"))         # returns the number of observations used 
                                                                            # in the comptation of the PCA (n for classic)
if(!isGeneric("getLoadings"))
    setGeneric("getLoadings", function(obj) standardGeneric("getLoadings")) 
if(!isGeneric("getEigenvalues"))
    setGeneric("getEigenvalues", function(obj) standardGeneric("getEigenvalues")) 
if(!isGeneric("getSdev"))
    setGeneric("getSdev", function(obj) standardGeneric("getSdev")) 
if(!isGeneric("getScores"))
    setGeneric("getScores", function(obj) standardGeneric("getScores")) 
if(!isGeneric("getPrcomp"))
    setGeneric("getPrcomp", function(obj) standardGeneric("getPrcomp"))     # return a prcomp() compatible object to use the 
                                                                            # available standard plots (i.e. S3 compatible obj.
                                                                            # with sdev, scale, scores, rotation
