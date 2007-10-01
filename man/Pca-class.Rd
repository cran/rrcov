\name{Pca-class}
\docType{class}
\alias{Pca-class}
\alias{getCenter,Pca-method}
\alias{getEigenvalues,Pca-method}
\alias{getLoadings,Pca-method}
\alias{getPrcomp,Pca-method}
\alias{getScores,Pca-method}
\alias{getSdev,Pca-method}
\alias{plot,Pca-method}
\alias{print,Pca-method}
\alias{show,Pca-method}

\title{Class "Pca" - virtual base class for all classic and robust PCA classes 
}
\description{  The class \code{Pca} searves as a base class for deriving all other 
classes representing the results of the classical and robust Principal 
Component Analisys methods
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{center}:}{Object of class \code{"vector"} the center of the data }
    \item{\code{loadings}:}{Object of class \code{"matrix"} the matrix 
    of variable loadings (i.e., a matrix whose columns contain the eigenvectors) }
    \item{\code{eigenvalues}:}{Object of class \code{"vector"} the eigenvalues }
    \item{\code{scores}:}{Object of class \code{"matrix"} the scores - the value 
    of the rotated data (the centred (and scaled if requested) data multiplied 
    by the \code{rotation} matrix) is returned.  Hence, \code{cov(scores)} 
    is the diagonal matrix \code{diag(eigenvalues)} }
    \item{\code{k}:}{Object of class \code{"numeric"} number of (choosen) principal components }
    \item{\code{sd}:}{Object of class \code{"Uvector"} Score distances within the robust PCA subspace }
    \item{\code{od}:}{Object of class \code{"Uvector"} Orthogonal distances to the robust PCA subspace }
    \item{\code{cutoff.sd}:}{Object of class \code{"numeric"} Cutoff value for the score distances}
    \item{\code{cutoff.od}:}{Object of class \code{"numeric"} Cutoff values for the orthogonal distances }
    \item{\code{flag}:}{Object of class \code{"Uvector"} The observations whose score distance is larger 
    than cutoff.sd or whose orthogonal distance is larger than cutoff.od can be considered 
    as outliers and receive a flag equal to zero.
    The regular observations receive a flag 1 }
    \item{\code{n.obs}:}{Object of class \code{"numeric"} the number of observations }
  }
}
\section{Methods}{
  \describe{
    \item{getCenter}{\code{signature(obj = "Pca")}: center of the data }
    \item{getEigenvalues}{\code{signature(obj = "Pca")}: the eigenvalues of the 
        covariance/correlation matrix, though the calculation is actually done 
        with the singular values of the data matrix) }
    \item{getLoadings}{\code{signature(obj = "Pca")}: returns the matrix of variable 
        loadings (i.e., a matrix whose columns contain the eigenvectors). 
        The function prcomp returns this matrix in the element rotation. }
    \item{getPrcomp}{\code{signature(obj = "Pca")}: returns an S3 object \code{prcomp} 
        for compatibility with the functions prcomp() and princomp(). Thus the 
        standard plots screeplot() and biplot() can be used}
    \item{getScores}{\code{signature(obj = "Pca")}:  returns the rotated data (the centred 
        (and scaled if requested) data multiplied by the loadings matrix). }
    \item{getSdev}{\code{signature(obj = "Pca")}: returns the standard deviations of the 
        principal components (i.e., the square roots of the eigenvalues of the 
        covariance/correlation matrix, though the calculation is actually done 
        with the singular values of the data matrix) }
    \item{plot}{\code{signature(x = "Pca")}: produces a distance plot (if \code{k=rank}) or 
        distance-distance plot (if\code{k<rank}) }
    \item{print}{\code{signature(x = "Pca")}: prints the results. The difference to the show()
        method is that additional parametesr ar epossible.}
    \item{show}{\code{signature(object = "Pca")}: prints the results }
     }
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{PcaClassic}}, \code{\link{PcaClassic-class}}, \code{\link{PcaRobust-class}}
}
\examples{
showClass("Pca")
}
\keyword{robust}
\keyword{multivariate}
