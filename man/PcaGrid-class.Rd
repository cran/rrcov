\name{PcaGrid-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaGrid-class}
\alias{getQuan,PcaGrid-method}

\title{Class "Pcaproj" - Robust PCA using PP - GRID search Algorithm  }
\description{
    Holds the results of an approximation of the PP-estimators for PCA using the grid search algorithm in the plane.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaProj", ...)} but the
    usual way of creating \code{PcaProj} objects is a call to the function
    \code{PcaProj()} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{call}, \code{center}, \code{loadings},
            \code{eigenvalues}, \code{scores}, \code{k},
            \code{sd}, \code{od}, \code{cutoff.sd}, \code{cutoff.od},
            \code{flag}, \code{n.obs}:}{
            from the \code{"\linkS4class{Pca}"} class.
        }
    }
}

\section{Extends}{
    Class \code{"\linkS4class{PcaRobust}"}, directly.
    Class \code{"\linkS4class{Pca}"}, by class \code{"\linkS4class{PcaRobust}"}, distance 2.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "PcaGrid")}: ... }
     }
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{PcaRobust-class}}, \code{\link{Pca-class}}, \code{\link{PcaClassic}}, \code{\link{PcaClassic-class}}
}

\examples{
showClass("PcaGrid")
}
\keyword{robust}
\keyword{multivariate}
