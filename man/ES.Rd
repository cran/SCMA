\name{ES}
\alias{ES}
\title{measure of effect size}
\description{Calculates the specified effect size measure.}
\usage{
ES(design, ES, data = read.table(file.choose(new = FALSE)))
}
\arguments{
  \item{design}{
Type of single-case design: "AB", "ABA", "ABAB", "CRD"(completely randomized design), "RBD" (randomized block design), "ATD" (alternating treatments design), "MBD" (multiple-baseline AB design) or "Custom" (user specified design).
}
  \item{ES}{
the type of effect size that has to be calculated: "SMD" (standardized mean difference), "SMDpool" (pooled standardized mean difference), "PND+" / "PND-" (percentage of nonoverlapping data, depending on the expected direction of the treatment effect), or "PEM+" / "PEM-" (percentage of data points exceeding the median, depending on the expected direction of the treatment effect)
}
  \item{data}{
File in which the data can be found. Default: a window pops up in which the appropriate file can be selected.
}
}
\details{
When using the default 'data' argument, a window will pop up to ask in what file the data can be found. This text file containing the data should consist of two columns for single-case phase and alternation designs: the first with the condition labels and the second with the obtained scores.
For multiple-baseline designs it should consist of these two columns for EACH unit. This way, each row represents one measurement occasion. It is important not to label the rows or columns.
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
data(AB)
ES(design = "AB", ES = "SMD", data = AB)}
\keyword{single-case designs}
\keyword{meta-analysis}
