\name{ES}

\alias{ES}

\title{Effect size measures for single-case data}

\description{The specified effect size measure is calculated: (pooled) standardized mean difference, percentage of nonoverlapping data (PND), or percentage of data points exceeding the median (PEM).}

\usage{ES(design, ES, data = read.table(file.choose(new = FALSE)))}

\arguments{
  \item{design}{Type of single-case design: \code{"AB"}, \code{"ABA"}, \code{"ABAB"}, \code{"CRD"}(completely randomized design), \code{"RBD"} (randomized block design), \code{"ATD"} (alternating treatments design), or \code{"MBD"} (multiple-baseline AB design).}

  \item{ES}{Effect size measure: the standardized mean difference where the standard deviation of the A phase is used as denominator (\code{"SMD"}), the pooled standardized mean difference which uses the square root of the pooled variance as denominator (\code{"SMDpool"}), the percentage of nonoverlapping data (\code{"PND+"} or \code{"PND-"}, depending on whether an increase or a decrease of behavior is expected following the intervention), or the percentage of data points exceeding the median (\code{"PEM+"} or \code{"PEM-"}).}

  \item{data}{File in which the data can be found. Default: a window pops up in which the file can be selected.}
}

\details{
When using the default \code{data} argument, a window will pop up to ask in what file the data can be found. This text file containing the data should consist of two columns for single-case phase and alternation designs: the first with the condition labels ("A" and "B" for studies with two conditions/phases and "A1", "A2", "B1" and "B2" when there are more phases/conditions) and the second with the obtained scores.
For multiple-baseline designs it should consist of these two columns for EACH unit. This way, each row represents one measurement occasion. It is important not to label the rows or columns.
}

\references{
Bulte, I., Van Den Noortgate, W., Onghena, P. (2010). An R package for the nonparametric meta-analysis of small-n educational studies. Annual Meeting of the American Educational Research Association. Denver, Colorado, April 30 - May 4.

\url{http://ppw.kuleuven.be/english/research/mesrg}
}

\examples{
data(AB)
ES(design = "AB", ES="SMD", data = AB)
}

\keyword{Single-case design}
\keyword{Effect size measure}
\keyword{Standardized mean difference}
\keyword{PND}
\keyword{PEM}