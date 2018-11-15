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
        
      \item  Added data sets \code{octane} and \code{diabetes}; data set \code{olitos} from package rrcovHD moved here.
      
      \item The 'pairs' plot is now available for classical covariance matrix
    }
  }
}
