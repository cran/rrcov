\name{lmom32}
\alias{lmom32}
\docType{data}
\title{ Hosking and Wallis Data Set, Table 3.2}
\description{
  The data on annual maximum streamflow at 18 sites with smallest
  drainage area basin in southeastern USA contains the sample
  L-moments ratios (L-CV, L-skewness and L-kurtosis) as
  used by Hosking and Wallis (1997) to illustrate the discordancy measure
  in regional freqency analysis (RFA).
}
\usage{data(lmom32)}

\format{
  A data frame with 18 observations on the following 3 variables.

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
    Approach Based on L-moments.} Cambridge University Press, p.49, Table 3.2
}
\references{
    Neykov, N.M., Neytchev, P.N., Van Gelder, P.H.A.J.M. and Todorov V. (2007),
    Robust detection of discordant sites in regional frequency analysis,
    \emph{Water Resources Research},
    to appear
}
\examples{
    data(lmom32)

    # plot a matrix of scatterplots
    pairs(lmom32,
          main="Hosking and Wallis Data Set, Table 3.3",
          pch=21,
          bg=c("red", "green3", "blue"))

    mcd<-CovMcd(lmom32)
    mcd
    plot(mcd, which="dist", class=TRUE)
    plot(mcd, which="dd", class=TRUE)
    sqrt(getDistance(mcd))
    cc <- Cov(lmom32)
    sqrt(getDistance(cc))
}
\keyword{datasets}
