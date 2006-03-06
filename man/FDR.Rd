\name{FDR}
\alias{FDR}
\alias{CDF}
\alias{CDFmix}
\alias{FDR.paired}
\alias{CDF.paired}
\alias{CDFmix.paired}
\title{Compute FDR for general scenarios}
\description{
\code{FDR} computes the false discovery rate for comparing gene expression 
between two groups of subjects when the distribution of the test statistic 
under the null and alternative hypothesis are both mixtures of t-distributions.
\code{CDF} and \code{CDFmix} calculate these mixtures.
}
\usage{
FDR(x, n1, n2, pmix, D0, p0, D1, p1, sigma)

CDF(x, n1, n2, D, p, sigma)
CDFmix(x, n1, n2, pmix, D0, p0, D1, p1, sigma)

FDR.paired(x, n, pmix, D0, p0, D1, p1, sigma)

CDF.paired(x, n, D, p, sigma)
CDFmix.paired(x, n, pmix, D0, p0, D1, p1, sigma)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of quantiles (two-sample t-statistics)}
  \item{n, n1, n2}{vector of sample sizes (as subjects per group)}
  \item{pmix}{the proportion of non-differentially expressed genes}
  \item{D0}{vector of effect sizes for the null distribution}
  \item{p0}{vector of mixing proportions for \code{D0}; must be the 
            same length as \code{D0} and sum to one}
  \item{D1}{vector of effect sizes for the alternative distribution}
  \item{p1}{vector of mixing proportions for \code{D1}, same as 
            \code{p0}}
  \item{D, p}{generic vectors of effect sizes and mixing proportions as above}            
  \item{sigma}{the standard deviation}
}
\details{
These functions are designed for a simple experimental setup, where we wish to
compare gene expression between two groups of subjects of size \code{n1} and
\code{n2} for an unspecified number of genes, using an equal-variance 
t-statistic. 

100\code{pmix}\% of the genes are assumed to be not differentially
expressed. The corresponding t-statistics follow a mixture of t-distributions;
this is more general than the usual central t-distribution, because we may want
to include genes with biologically small effects under the null hypothesis
(Pawitan et al., 2005). The other 100(1-\code{pmix})\% genes are assumed to be differentially expressed; their t-statistics are also mixtures of t-distributions. 

The mixture proportions of t-distributions under the null and alternative
hypothesis are specified via \code{p0} and \code{p1}, respectively. The
individual t-distributions are specified via the means \code{D0} and \code{D1} and the standard deviation \code{sigma} of the underlying data (instead of the mathematically more obvious, but less intuitive non centrality parameters). As the underlying data are the logarithmized expression values, \code{D0} and \code{D1} can be interpreted as average log-fold change between conditions, measured in units of \code{sigma}. See Examples.

\code{CDF} computes the cumulative distribution function for a mixture of
t-distributions based on means \code{D} and standard deviation \code{sigma} with
mixture proportions \code{p}. This function is the work horse for \code{CDFmix}.

Note that the base functions (\code{FDR}, \code{CDFmix}, \code{CDF}) assume two groups of experimental units; the \code{.paired} functions provide the same functionality for one group of paired observations. 

The distribution functions call \code{pt} for computation; correspondingly, the quantiles \code{x} and all arguments that define degrees of freedom and non centrality parameters (\code{n1}, \code{n2}, \code{D0}, \code{D1}, \code{sigma}) can be vectors, and will be recycled as necessary.
}
\value{
The appropriate vector of FDRs or probabilities.
}
\references{
Pawitan Y, Michiels S, Koscielny S, Gusnanto A, Ploner A. (2005) False Discovery Rate, Sensitivity and Sample Size for Microarray Studies. \emph{Bioinformatics}, 21, 3017-3024.
}
\author{Y. Pawitan and A. Ploner}
%\note{ ~~further notes~~ }
\seealso{\code{\link{TOC}}, \code{\link{samplesize}}}
\examples{
# FDR for H0: 'log fold change is zero'
#     vs. H1: 'log fold change is -1 or 1' 
#             (ie two-fold up- or down regulation) 
FDR(1:6, n1=10, n2=10, pmix=0.90, D0=0, p0=1, 
    D1=c(-1,1), p1=c(0.5, 0.5), sigma=1)

# Include small log fold changes in the H0
# Naturally, this increases the FDR
FDR(1:6, n1=10, n2=10, pmix=0.90, D0=c(-0.25,0, 0.25), p0=c(1/3,1/3,1/3), 
    D1=c(-1,1), p1=c(0.5, 0.5), sigma=1)

# Consider an asymmetric alternative
# 10 percent of the regulated genes are assumed to be four-fold upregulated
FDR(1:6, n1=10, n2=10, pmix=0.90, D0=0, p0=1, 
    D1=c(-1,1,2), p1=c(0.45, 0.45, 0.1), sigma=1)
}
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{distribution}% __ONLY ONE__ keyword per line
