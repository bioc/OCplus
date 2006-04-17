\name{EOC}
\alias{EOC}
\alias{FDRp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimated or empirical FDR, sensitivity, etc as a function of cutoff level}
\description{
\code{EOC} computes and optionally plots the estimated operating characteristics for data from a microarray experiment with two groups of subjects. The false discovery rate (FDR) is estimated based on random permutations of the data and plotted against the cutoff level on the t-statistic; a curve for the classical sensitivity can be added. Different curves for different proportions of non-differentially expressed genes can be compared in the same plot, and the sample size per group can be varied between plots.

\code{FDRp} is the function that does the underlying hard work and requires package \code{multtest}.
}
\usage{
EOC(xdat, grp, p0, paired = FALSE, nperm = 25, seed = NULL, plot = TRUE, ...)

FDRp(xdat, grp, test = "t.equalvar", p0, nperm, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdat}{the matrix of expression values, with genes as rows and samples as columns} 
  \item{grp}{a grouping variable giving the class membership of each sample, i.e. each column in \code{xdat}; for \code{EOC}, this can be any type of variable, as long as it has exactly two distinct values, whereas \code{FDRp} expects to see only 0s and 1s, see Details.}
  \item{p0}{if supplied, an estimate for the proportion of non-differentially expressed genes; if not supplied, the routine will estimate it, see Details.}
  \item{paired}{logical value indicating whether this is independent sample
  situation (default) or a paired sample situation. Note that paired samples need to follow each other in the data matrix (as in 010101...} when \code{paired=TRUE}.
  \item{nperm}{number of permutations for establishing the null distribution of the t-statistic}
  \item{test}{the type of test to use, see \code{mt.teststat}; when called from \code{EOC}, this is always the default.}  
  \item{seed}{the random seed from which the permutations are started}
  \item{plot}{logical value indicating whether to do the plot}
  \item{\dots}{graphical parameters, passed to \code{plot.FDR.result}}
}
\details{
\code{EOC} is the empirical counterpart of the function \code{TOC}. It estimates the FDR and sensitivity for a given data set of expression values measured on subjects in two groups. The FDR is estimated locally based on the empirical Bayes approach outlined by Efron et al., see References. \code{FDRp} implements the details of this method; this requires among other things the permutation distribution of the t-statistic, which is calculated via a call to function \code{mt.teststat} of package \code{multtest}. This explains why both functions barf at missing values in the expression data.

Note that \code{p0} is by default estimated from the data, as originally suggested by Efron et al. so as to make ratio between the densities of the observed distribution of t-statistics and the permutation distribution smaller than 1; alternatively, the user can supply his own guesstimate of the proportion of non-differentially expressed genes in the data.

Note also that \code{FDRp} keeps all permuations in the memory during compuations. For a large number of genes, this will limit the number of possible permuations. 
}
\value{
For \code{EOC}, an object of class \code{FDR.result}, which inherits from class \code{data.frame}. The three columns list for each gene its t-statistic, the estimated FDR (two-sided), and the estimated sensitivity. Additionally, the object carries an attribute \code{param}, which is a list with four entries: \code{p0}, the assumed proportion of non-differentially expressed genes used in calculating the FDR; \code{p0.est}, a logical value indicating whether \code{p0} was estimated or user-supplied; \code{statistic} indicates how the t-statistic was computed, i.e. how its sign should be interpreted in terms of relative over- or under expression, and a logical flag \code{paired} to indicate whether a paired t-statistic was used. 

\code{FDRp} returns a list with essentially the same elements, plus additionally the values of the observed and permuted distribution of the t-statistics for each gene.
}
\references{
Pawitan Y, Michiels S, Koscielny S, Gusnanto A, Ploner A (2005) False Discovery Rate, Sensitivity and Sample Size for Microarray Studies. \emph{Bioinformatics}, 21, 3017-3024.

Efron B, Tibshirani R, Storey JD, Tusher V. (2001) Empirical Bayes Analysis of a Microarray Experiment. \emph{JASA}, 96(456), p. 1151-60.
}
\author{Y. Pawitan and A. Ploner}
\note{Both the curve labels and the legend may be squashed if the plotting device is too small. Increasing the size of the device and re-plotting should improve readability.}
\seealso{\code{\link{plot.FDR.result}}, \code{\link{OCshow}}, \code{\link[multtest]{mt.teststat}}}
\examples{
# We simulate a small example with 5 percent regulated genes and
# a rather large effect size
set.seed(2003)
xdat = matrix(rnorm(50000), nrow=1000)
xdat[1:25, 1:25] = xdat[1:25, 1:25] - 2
xdat[26:50, 1:25] = xdat[26:50, 1:25] + 2
grp = rep(c("Sample A","Sample B"), c(25,25))

# The default, with legend
ret = EOC(xdat, grp, legend=TRUE)
# Look at the results: yes
ret[1:10,]
which(ret$FDR<0.05)
# Extra information
attr(ret,"param")

# Run the same data with different permutations: fairly stable, but with
# different p0
ret = EOC(xdat, grp, seed=2000)
which(ret$FDR<0.07)

# Misspecify the p0: not too bad here
ret = EOC(xdat, grp, p0=0.99)
which(ret$FDR<0.01)

# We simulate data in a paired setting
# Note the arrangement of the columns
set.seed(2004)
xdat = matrix(rnorm(50000), nrow=1000)
ndx1 = seq(1,50, by=2)
xdat[1:25, ndx1] = xdat[1:25, ndx1] - 2
xdat[26:50, ndx1] = xdat[26:50, ndx1] + 2
grp = rep(c("Sample A","Sample B"), 25)

ret = EOC(xdat, grp, paired=TRUE)
which(ret$FDR<0.05)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
