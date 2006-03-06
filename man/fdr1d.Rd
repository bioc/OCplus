\name{fdr1d}
\alias{fdr1d}
\title{Compute classical local false discovery rate}
\description{
Calculates the classical local false discovery rate for multiple 
parallel t-statistics. 
}
\usage{
fdr1d(xdat, grp, test, p0, nperm = 100, nr = 50, seed = NULL, null = NULL, 
      zlim = 1, sv2 = 0.01, err = 1e-04, verb = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdat}{the matrix of expression values, with genes as rows and samples as columns}
  \item{grp}{a grouping variable giving the class membership of each sample, i.e. each column in \code{xdat}}
  \item{test}{a function that takes \code{xdat} and \code{grp} as the first two
  arguments and returns the test statistic; by default, two-sample t-statistics are calculated.}
  \item{p0}{if supplied, an estimate for the proportion of non-differentially expressed genes; if not supplied, the routine will estimate it, see Details.}
  \item{nperm}{number of permutations for establishing the null distribution of the t-statistic}
  \item{nr}{the number of equidistant breaks into which the range of test statistics is broken for calculating the fdr.}
  \item{seed}{if specified, the random seed from which the permuations are
  started}
  \item{null}{optional argument for passing in a pre-calculated null distribution, see Details.}
  \item{zlim}{if no \code{p0} is specified, the ratio of densities in the range of test statistics between \code{-zlim} and \code{zlim} will be used to estimate the proportion of non-differentially expressed genes; ignored if \code{p0} is specified.}
  \item{sv2}{positive number controlling the initial degree of smoothing for the densities involved, with smaller values indicating more smoothing; see Details.}
  \item{err}{positive number controlling the convergence of the smoothing procedure, with smaller values implying more iterations; see Details. }
  \item{verb}{logical value indicating whether provide extra information.}
  \item{\dots}{extra arguments to function \code{test}.}
}
\details{
This function calculates the local false discovery rate (fdr, as opposed to global FDR) introduced by Efron et al., 2001. The underlying model assumes that for a given grouping of samples, we study a mixture of differentially expressed (DE) and non-DE genes, and that consequently, the observed distribution of test statistics is a mixture of test statistics under the alternative and the null statistic, respectively. The densities involved are estimated nonparametrically and smoothed, using a permutation argument for the null distribution.

By default, the null distribution is generated and stored only within the function, and is not available outside. If someone wants to study the null distribution, or wants to re-use the same null distribution efficiently while e.g. varying the smoothing parameter, the argument \code{null} allows the use of an externally generated null distribution, created e.g. using the \code{PermNull} function.

Theoretically, the function should support any kind of function along the lines of \code{tstatistics}, however, this has not been tested, and the current setup is very much geared towards t-tests.

We use non-Gaussian mixed model smoothing for Possion counts for smoothing the density estimates, see Pawitan, 2001, and \code{smooth1d}.
}
\value{
Basically, a data frame with one row per gene and two columns: \code{tstat}, the test statistic, and \code{fdr.local}, the local false discovery rate. This data frame has the additional class attributes \code{fdr1d.result} and \code{fdr.result}, see Examples. This is the bad old S3 class mechanism employed to provide plot and summary functions. 

Additional information is provided by a \code{param} attribute, which is a list with the following entries:
  \item{p0}{the proportion of non-differentially expressed genes used when calculating the fdr.}
  \item{p0.est}{a logical value indicating whether \code{p0} was estimated from the data or supplied by the user.}
  \item{fdr}{the smoothed fdr values calculated for the original intervals.}
  \item{xbreaks}{vector of breaks for the test statistic defining the interval for calculation.}
}
\references{
Efron B, Tibshirani R, Storey JD, Tusher V (2001) Empirical Bayes Analysis of a Microarray Experiment. \emph{JASA}, 96(456), p. 1151-60.

Pawitan Y.(2001) \emph{In All Likelihood}, Oxford University Press, ch. 18.11
}
\author{A. Ploner}
\seealso{\code{\link{plot.fdr1d.result}}, \code{\link{summary.fdr.result}}, \code{\link{OCshow}}, \code{\link{tstatistics}}, \code{\link{smooth1d}}, \code{\link{fdr2d}}, \code{\link{PermNull}}}
\examples{
# We simulate a small example with 5 percent regulated genes and
# a rather large effect size
set.seed(2000)
xdat = matrix(rnorm(50000), nrow=1000)
xdat[1:25, 1:25] = xdat[1:25, 1:25] - 1
xdat[26:50, 1:25] = xdat[26:50, 1:25] + 1
grp = rep(c("Sample A","Sample B"), c(25,25))

# A default run
res1d = fdr1d(xdat, grp)
res1d[1:20,]

# Looking at the results
summary(res1d)
plot(res1d)
res1d[res1d$fdr<0.05, ]

# Averaging estimates the global FDR for a set of genes
ndx = abs(res1d$tstat) > 3
mean(res1d$fdr[ndx])

# Extra information
class(res1d)
attr(res1d,"param")

}
\keyword{htest}% at least one, from doc/KEYWORDS

