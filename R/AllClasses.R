setClass("PsiFun", representation(
                               n = "numeric",
                               p = "numeric",
                               r = "numeric",
                               alpha = "numeric",
                               c1 = "numeric"))

setClass("PsiBwt", representation(
                                M = "numeric"), 
                   contains="PsiFun")
                   
