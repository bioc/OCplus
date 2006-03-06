\name{average.fdr}
\alias{average.fdr}
\title{Average a two-dimensional local false discovery rate}
\description{
This function averages two-dimensional local false discovery rates as computed by \code{fdr2d} for binned values of the first test statistic and across the values of the second test statistic. The result can easily be plotted and should be comparable to the one-dimensional fdr as provided by \code{fdr1d}, provided that the smoothing parameters were chosen suitably.
}
\usage{
average.fdr(x, breaks)
}
\arguments{
  \item{x}{an object returned by \code{fdr2d}.}
  \item{breaks}{breaks defining intervals into which the first test statistic is binned; by default the same values that were used by \code{fdr2d}.}
}
\details{
Assuming that we have smoothed the estimate suitably and have chosen the proportion of non-dffierentially expressed genes suitably, we should get very much the same results from \code{fdr2d} as from \code{fdr1d} when we average across the logarithmized standard errors, see Examples.

The averaging is done across the estimated values for the actual genes; this corresponds to a weighted mean of the smoothed estimates on a grid, where the weight is proportional to cell frequencies. 

Note that it is usuually easier to get a good match in the tails of the curves than in the center, which is okay, as this is where we want to estimate the fdr reliably.
}
\value{
A matrix with two columns \code{tstat} and \code{fdr.local}. 
}
\references{
Ploner A, Calza S, Gusnanto A, Pawitan Y (2005) Multidimensional local false discovery rate for micorarray studies. \emph{Submitted Manuscript}.
}
\author{A. Ploner}
\seealso{\code{\link{fdr2d}}, \code{\link{fdr1d}}}
\examples{
# Create res1d
example(fdr1d)

# Compute fdr2d using the p0
res2d = fdr2d(xdat, grp, p0=p0(res1d))

# Show it
par(mfrow=c(2,1))
plot(res1d, main="fdr1d and averaged fdr2d")
lines(average.fdr(res2d), col="red")
plot(res2d, grid=TRUE, main="fdr2d is averaged across columns")
}
\keyword{utilities}% at least one, from doc/KEYWORDS

