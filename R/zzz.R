.First.lib <- function(lib, pkg) {

    where <- match(paste("package:", pkg, sep = ""), search())
    ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    ver <- as.character(ver)
    title <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Title")
    title <- as.character(title)
    cat(paste(title, " (version ", ver, ")\n", sep = ""))
    
    library.dynam("rrcov", pkg, lib)
}
