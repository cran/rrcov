% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datadoc.R
\docType{data}
\name{ionosphere}
\alias{ionosphere}
\title{Johns Hopkins University Ionosphere database.}
\format{
A data frame with 351 rows and 33 variables: 32 measurements and one
 (the last, \code{Class}) grouping variable: 225 \code{'good'} and 126 \code{'bad'}.

 The original dataset at UCI contains 351 rows and 35 columns. The first 34
 columns are features, the last column contains the classification label of
 'g' and 'b'. The first feature is binary and the second one is only 0s,
 one grouping variable - factor with labels 'good' and 'bad'.
}
\source{
Source: Space Physics Group; Applied Physics Laboratory; Johns Hopkins University; Johns Hopkins Road; Laurel; MD 20723

 Donor: Vince Sigillito (vgs@aplcen.apl.jhu.edu)

 The data have been taken from the UCI Repository Of Machine Learning Databases at
 \url{https://archive.ics.uci.edu/ml/datasets/ionosphere}

 This data set, with the original 34 features is available in the package \pkg{mlbench}
 and a different data set (refering to the same UCI repository) is available in
 the package \code{dprep} (archived on CRAN).
}
\usage{
data(ionosphere)
}
\description{
''This radar data was collected by a system in Goose Bay, Labrador.  This
  system consists of a phased array of 16 high-frequency antennas with a
  total transmitted power on the order of 6.4 kilowatts.  The targets
  were free electrons in the ionosphere.
  "good" radar returns are those showing evidence of some type of structure
  in the ionosphere.  "bad" returns are those that do not; their signals pass
  through the ionosphere.
  Received signals were processed using an autocorrelation function whose
  arguments are the time of a pulse and the pulse number.  There were 17
  described by 2 attributes per pulse number, corresponding to the complex
  values returned by the function resulting from the complex electromagnetic
  signal.'' [UCI archive]
}
\examples{
 data(ionosphere)
 pairs(ionosphere[, 1:6])
}
\references{
Sigillito, V. G., Wing, S. P., Hutton, L. V., and Baker, K. B. (1989).
     Classification of radar returns from the ionosphere using neural
     networks. Johns Hopkins APL Technical Digest, 10, 262-266.
}
