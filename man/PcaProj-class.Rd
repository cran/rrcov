\name{PcaProj-class}
\docType{class}
\alias{PcaProj-class}
\alias{getQuan,PcaProj-method}

\title{Class "Pcaproj" - Robust PCA using PP - Croux and Ruiz-Gazen (2005) algorithm }
\description{
    Holds the results of an approximation of the PP-estimators for PCA by a fast and simple algorithm: Croux and Ruiz-Gazen (2005) algorithm.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaProj", ...)} but the 
    usual way of creating \code{PcaProj} objects is a call to the function
    \code{PcaProj()} which serves as a constructor.
}
\section{Slots}{
     \describe{
    \item{\code{delta}:}{Accuracy parameter}
    \item{\code{quan}:}{Object of class \code{"numeric"} The quantile h used throughout the algorithm }
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{center}:}{Object of class \code{"vector"} the center of the data }
    \item{\code{loadings}:}{Object of class \code{"matrix"} the matrix 
    of variable loadings (i.e., a matrix whose columns contain the eigenvectors) }
    \item{\code{eigenvalues}:}{Object of class \code{"vector"} the eigenvalues }
    \item{\code{scores}:}{Object of class \code{"matrix"} the scores - the value 
    of the projected on the space of the principal components data (the centred 
    (and scaled if requested) data multiplied 
    by the \code{loadings} matrix) is returned.  Hence, \code{cov(scores)} 
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

\section{Extends}{
Class \code{"\linkS4class{PcaRobust}"}, directly.
Class \code{"\linkS4class{Pca}"}, by class "PcaRobust", distance 2.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "PcaProj")}: ... }
     }
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{PcaRobust-class}}, \code{\link{Pca-class}}, \code{\link{PcaClassic}}, \code{\link{PcaClassic-class}}
}

\examples{
showClass("PcaProj")
}
\keyword{robust}
\keyword{multivariate}
