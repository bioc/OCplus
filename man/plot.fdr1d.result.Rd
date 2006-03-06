\name{plot.fdr1d.result}
\alias{plot.fdr1d.result}
\title{Plot univariate local false discovery output}
\description{
A plotting function for \code{fdr1d}.
}
\usage{
plot.fdr1d.result(x, add = FALSE, grid = FALSE, rug = TRUE, 
                  xlab = "t-Statistic", ylab = "fdr", lcol = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{output from \code{fdr1d}}
  \item{add}{logical value indicating whether to create a new plot or add to an existing one}
  \item{grid}{logical value indicating whether to show the intervals used for calculating the fdr.}
  \item{rug}{logical value indicating whether to add a 1D scatterplot showing the observed test statistics}
  \item{xlab, ylab}{the usual axis labels}
  \item{lcol}{the color of the lines}
  \item{\dots}{extra options passed to \code{plot.default}.}
}
\author{A Ploner}
\seealso{\code{\link{fdr1d}}}
\examples{
example(fdr1d)
plot(res1d, grid=TRUE, rug=FALSE)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{aplot}% __ONLY ONE__ keyword per line
