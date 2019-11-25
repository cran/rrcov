% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("C:/projects/statproj/R/rrcov/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{rrcov}}
\encoding{UTF-8}

\section{CHANGES in rrcov VERSION 1.4-7 (2018-11-14)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added Minimum Regularized Covariance Determinant estimator (MRCD) 
        (Boudt et al. 2018) - 
        functions \code{CovMrcd()} and \code{CovControlMrcd()}
        
      \item  Added data set \code{octane}; data set \code{olitos} from package rrcovHD moved here.
      
      \item The 'pairs' plot is now available for classical covariance matrix
    }
  }
}
\section{CHANGES in rrcov VERSION 1.4-8 (2019-11-15)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item  Added data set \code{diabetes}; data set \code{soil} from package \code{rrcovHD} moved here.
      \item  Linear and quadratic discriminant analysis can use the MRCD estimates.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue with CovControlMcd(nsamp="deterministic") - this would not work, 
        because nsamp was defined in the class definition as "numeric". Now it is "Cnumeric" - union of "character" and "numeric'.
      \item Corrected the code for Minimum Regularized Covariance Determinant estimator (MRCD) - 
        \code{CovMrcd()} - the step of adjusting the eignevalues in r6pack() is excluded 
        now because it has no effect when n > p.
    }
  }
}
