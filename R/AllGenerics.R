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

setGeneric("getCenter", function(obj) standardGeneric("getCenter")) 
setGeneric("getCov", function(obj) standardGeneric("getCov")) 
setGeneric("getCorr", function(obj) standardGeneric("getCorr")) 
setGeneric("getData", function(obj) standardGeneric("getData")) 
setGeneric("getDistance", function(obj) standardGeneric("getDistance")) 
setGeneric("getEvals", function(obj) standardGeneric("getEvals")) 

setGeneric("estimate", function(obj, x, ...) standardGeneric("estimate")) 
