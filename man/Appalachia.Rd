\name{Appalachia}
\alias{Appalachia}
\docType{data}
\title{ Annual maximum streamflow in central Appalachia}
\description{
  The data on annual maximum streamflow at 104 gaging stations
  in the central Appalachia region of the United States contains
  the sample L-moments ratios (L-CV, L-skewness and L-kurtosis) as
  used by Hosking and Wallis (1997) to illustrate regional
  freqency analysis (RFA).
}
\usage{data(Appalachia)}

\format{
  A data frame with 104 observations on the following 3 variables.

\describe{
    \item{\code{L-CV}}{L-coefficient of variation}
    \item{\code{L-skewness}}{L-coefficient of skewness}
    \item{\code{L-kurtosis}}{L-coefficient of kurtosis}
  }

}

\details{
The sample L-moment ratios (L-CV, L-skewness and L-kurtosis) of a
site are regarded as a point in three dimensional space.
}
\source{
    Hosking, J. R. M. and J. R. Wallis (1997), \emph{Regional Frequency Analysis: An
    Approach Based on L-moments.} Cambridge University Press, p.175--185
}
\references{
    Neykov, N.M., Neytchev, P.N., Van Gelder, P.H.A.J.M. and Todorov V. (2007),
    Robust detection of discordant sites in regional frequency analysis,
    \emph{Water Resources Research},
    to appear
}
\examples{
    data(Appalachia)

    # plot a matrix of scatterplots
    pairs(Appalachia,
          main="Appalachia data set",
          pch=21,
          bg=c("red", "green3", "blue"))

    mcd<-CovMcd(Appalachia)
    mcd
    plot(mcd, which="dist", class=TRUE)
    plot(mcd, which="dd", class=TRUE)

}
\keyword{datasets}
