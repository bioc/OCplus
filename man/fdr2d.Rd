\name{fdr2d}
\alias{fdr2d}
\title{Compute two-dimensional local false discovery rate}
\description{
This function calculates the local false discovery rate for a two-sample problem using a bivariate test statistic, consisting of classical t-statistics and the corresponding logarithmized standard error.
}
\usage{
fdr2d(xdat, grp, test, p0, nperm = 100, nr = 15, seed = NULL, null = NULL,
	  constrain = TRUE, smooth = 0.2, verb = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdat}{the matrix of expression values, with genes as rows and samples as columns}
  \item{grp}{a grouping variable giving the class membership of each sample, i.e. each column in \code{xdat}}
  \item{test}{a function that takes \code{xdat} and \code{grp} as the first two
  arguments and returns the bivariate test statistics as two-column matrix; by default, two-sample t-statistics and logrithmized standard errors are calculated.}
  \item{p0}{if supplied, an estimate for the proportion of non-differentially expressed genes; if not supplied, the routine will estimate it, see Details.}
  \item{nperm}{number of permutations for establishing the null distribution of the t-statistic}
  \item{nr}{the number of equidistant breaks for the range of each test statistic; fdr values are calculated on the resulting (nr-1) x (nr-1) grid of cells.}
  \item{seed}{if specified, the random seed from which the permuations are
  started}
  \item{null}{optional argument for passing in a pre-calculated null distribution, see Examples.}
  \item{constrain}{logical value indicating whether the estimated fdr should be constrained to be monotonously decreasing with the absolute size of the t-statistic (more generally, the first test statistic).}
  \item{smooth}{a numerical value between 0.01 and 0.99, indicating which percentage of the available degrees of freedom are used for smoothing the fdr estimate; larger values indicate more smoothing.}
  \item{verb}{logical value indicating whether provide extra information.}
  \item{\dots}{extra arguments to function \code{test}.}
}
\details{
This routine computes a bivariate extension of the classical local false discovery rate as available through function \code{fdr1d}. Consequently, many arguments have identical or similar meaning. Specifically for \code{fdr2d}, \code{nr} specifies the number of equidistant breaks defining a two-dimensional grid of cells on which the bivariate test statistics are counted; argument \code{constrain} can be set to ensure that the estimated fdr is decreasing with increasing absolute value of the t-statistic; and argument \code{smooth} specifies the degree of smoothing when estimating the fdr.

Note that while \code{fdr2d} might be used for any suitable pair of test statistics, it has only been tested for the default pair, and the smoothing procedure specifically is optimized for this situation.

Note also that the estimation of the proportion \code{p0} directly from the data may be quite unstable and dependant on the degree of smoothing; too heavy smoothing may even lead to estimates greater than 1. It is usually more stable 	use an estimate of \code{p0} provided by \code{fdr1d}.

Note that \code{fdr1d} can also be used to check the degree of smoothing, see \code{average.fdr}.
}
\value{
Basically, a data frame with one row per gene and three columns: \code{tstat}, the test statistic, \code{logse}, the corresponding logarithmized standard error, and \code{fdr.local}, the local false discovery rate. This data frame has the additional class attributes \code{fdr2d.result} and \code{fdr.result}, see Examples. This is the bad old S3 class mechanism employed to provide plot and summary functions. 

Additional information is provided by a \code{param} attribute, which is a list with the following entries:
  \item{p0}{the proportion of non-differentially expressed genes used when calculating the fdr.}
  \item{p0.est}{a logical value indicating whether \code{p0} was estimated from the data or supplied by the user.}
  \item{fdr}{the matrix of smoothed fdr values calculated on the original grid.}
  \item{xbreaks}{vector of breaks for the first test statistic.}
  \item{ybreaks}{vector of breaks for the second test statistic.}  
}
\references{
Ploner A, Calza S, Gusnanto A, Pawitan Y (2005) Multidimensional local false discovery rate for micorarray studies. \emph{Submitted Manuscript}.
}
\author{A Ploner and Y Pawitan}
\seealso{\code{\link{plot.fdr2d.result}}, \code{\link{summary.fdr.result}}, \code{\link{OCshow}}, \code{\link{fdr1d}}, \code{\link{average.fdr}}}
\examples{
# We simulate a small example with 5 percent regulated genes and
# a rather large effect size
set.seed(2000)
xdat = matrix(rnorm(50000), nrow=1000)
xdat[1:25, 1:25] = xdat[1:25, 1:25] - 1
xdat[26:50, 1:25] = xdat[26:50, 1:25] + 1
grp = rep(c("Sample A","Sample B"), c(25,25))

# A default run
res2d = fdr2d(xdat, grp)
res2d[1:20,]

# Looking at the results
summary(res2d)
plot(res2d)
res2d[res2d$fdr<0.05, ]

# Extra information
class(res2d)
attr(res2d,"param")
}
\keyword{htest}% at least one, from doc/KEYWORDS

