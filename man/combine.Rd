\name{combine}
\alias{combine}
\title{Statistical combination of p-values}
\description{Calculates a general p-value by statistically combining the p-values of a number of independent studies, to determine whether a general significant result is obtained.}
\usage{
combine(method, pvalues = read.table(file.choose(new = FALSE)))
}
\arguments{
  \item{method}{
Indicates which combining function should be used: "x" (multiplicative) or "+" (additive)
}
  \item{pvalues}{
File in which the p-values can be found. Default: a window pops up in which the appropriate file can be selected.
}
}
\details{
When using the default 'pvalues' argument, a window will pop up to ask in what file the pvalues can be found. This text file containing the pvalues should consist of 1 column with all the obtained pvalues.
}
\references{
Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467-478.

Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477-485.

http://ppw.kuleuven.be/cmes/SCRT-R.html
}
\author{
Isis Bulte
}
\examples{
data(P)
combine(method="+",pvalues=P)}
\keyword{single-case designs}
\keyword{meta-analysis}
