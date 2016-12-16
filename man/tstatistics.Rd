\name{tstatistics}
\alias{tstatistics}
\alias{PermNull}
\title{Compute multiple parallel t-statistics}
\description{
\code{tstatistics} computes either two-sample or paired t-statistics for a bunch of variables measured on the same objects, e.g. genewise t-statistics for a microarray experiment. \code{PermNull} uses \code{tstatistics} to generate a permutation distribution.
}
\usage{
tstatistics(xdat, grp, logse = FALSE, paired = FALSE)

PermNull(xdat, grp, nperm = 100, seed = NULL, logse = FALSE, paired=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdat}{the matrix of expression values, with genes (or variables) as rows and samples as columns.}
  \item{grp}{a grouping variable giving the class membership of each sample, i.e. each column in \code{xdat}, see Details.}
  \item{nperm}{number of permutations for establishing the null distribution of the t-statistic}
  \item{seed}{random number generator seed for initializing the permutations from a known starting point.}
  \item{logse}{logical flag indicating whether to return the logarithmized standard errors, too.}
  \item{paired}{indicates whether to use two-sample or paired t-statistic.}
}
\details{
\command{tstatistics}  is a fairly fast replacement for function \command{mt.teststat} in package \pkg{multtest}, which is written exlusively in \R and does not require loading half the Bioconductor infrastructure packages before doing anything. As such, it is used for computing the default test statistics by \command{fdr1d} and \command{fdr2d}.

Note that for the paired test, \command{tstatistics} requires the same data structure as \command{mt.teststat}: columns belonging to the same pair must be consecutive (though not necessarily in the same order throughout, as \option{grp} will indicate the order). The function checks for this and barfs if it does not hold.

\command{PermNull} returns the t-statistics and optionally the logarithmized standard errors of the mean for a specified number of permutations.

Both functions are not especially economic in using memory, and collecting the whole set of permutations like \command{PermNull} does instead of binning and counting them directly as they come is inherently wasteful.
}
\value{
A data frame with first column \code{tstat} and optionally (if \code{logse=TRUE}) a second column \code{logse}. \code{tstat} returns the same number of test statistics as  rows in \code{xdat} and in the same order, \code{PermNull} does the same for consecutive permuations of the grouping variable \code{grp}. 

If the argument \code{seed} is specified, \code{PermNull} adds an attribute of the same name to the returned data frame. 
}
\author{A. Ploner}
\seealso{\code{\link{fdr1d}}, \code{\link{fdr2d}}, examples in \code{\link{tMixture}}}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{utilities} % __ONLY ONE__ keyword per line
