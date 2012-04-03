\name{combine}

\alias{combine}

\title{Probability combining: nonparametric combination of p-values}

\description{The \code{combine} function contains two methods for the nonparametric combination of p-value: the multiplicative approach using Pearson's formula, and the additive approach using Edgington's formula.}

\usage{combine(method, pvalues = read.table(file.choose(new = FALSE)))}

\arguments{
  \item{method}{Probability combining method: multiplicative approach (\code{method = "x"}) or additive approach (\code{method = "+"}).}

  \item{pvalues}{File in which the p-values to be combined can be found. Default: a window pops up in which the file can be selected.}
}

\details{
When using the default \code{pvalues} argument, a window will pop up to ask in what file the p-values can be found. This text file containing the p-values should consist of one column (all p-values below one another), without column or row labels.}

\references{
Bulte, I., Van Den Noortgate, W., Onghena, P. (2010). An R package for the nonparametric meta-analysis of small-n educational studies. Annual Meeting of the American Educational Research Association. Denver, Colorado, April 30 - May 4.

\url{http://ppw.kuleuven.be/english/research/mesrg}
}

\examples{
data(P)
combine(method = "+", pvalues = P)
}

\keyword{Probability combining}
\keyword{Meta-analysis}