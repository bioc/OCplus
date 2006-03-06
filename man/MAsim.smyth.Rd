\name{MAsim.smyth}
\alias{MAsim.smyth}
\alias{MAsim}
\alias{MAsim.var}
\alias{MAsim.real}
\title{Simulate two-sample microarray data}
\description{
These functions simulate two-sample microarray data from various different models.
}
\usage{
MAsim(ng = 10000, n = 10, n1 = n, n2 = n, D = 1, p0 = 0.9, sigma = 1)


MAsim.var(ng = 10000, n = 10, n1 = n, n2 = n, D = 1, p0 = 0.9) 

MAsim.smyth(ng = 10000, n = 10, n1 = n, n2 = n, p0 = 0.9, d0 = 4, 
            s2_0 = 4, v0 = 2)

MAsim.real(xdat, grp, n, n1, n2, D = 1, p0 = 0.9, replace = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ng}{number of genes}
  \item{n, n1, n2}{number of samples per group; by default balanced, except for \code{MAsim.real}.}
  \item{p0}{proportion of differentially expressed genes}
  \item{D}{effect size for differentially expressed genes, in units of the gene-specific standard deviation (\code{sigma} in \code{MAsim}).}
  \item{sigma}{standard deviation, constant for all genes}
  \item{d0, s2_0, v0}{prior parameters for effect size and variability across genes in Smyth's model, see Details.}
  \item{xdat, grp}{expression data and grouping variable for an existing microarray data set, as specified in \code{EOC}.}
  \item{replace}{logical switch indicating whether to sub-sample (\code{replace=FALSE}) or bootstrap (\code{replace=TRUE}) from the existing data. Note that the specified group-sizes have to be smaller than the real group sizes in case of sub-sampling.}
}
\details{
\code{MAsim} simulates normal data with constant standard deviation \code{sigma} across genes and fixed effect size \code{D}; the sign of the effect is equally and randomly split between up- and down-regulation, and effects are added to the second group. \code{MAsim.var} does the same, but instead of relying on a fixed variance across genes, it simulates gene-specific variances from a standard exponential distribution.

\code{MAsim.smyth} simulates from the model suggested in Smyth (2004), using a normal error distribution. The variances are assumed to follow an inverse chisquared distribution with \code{d0} degrees of freedom and are scaled by \code{s2_0}; consequently, large values of \code{d0} lead to similar gene-wise variances across genes, whereas small values lead to very different variances between genes. The effect sizes for differentially expressed genes are assumed to follow a normal distribution with mean zero and variance \code{v0} times the previously simulated gene-specific variance; consequently, large values of \code{v0} lead to large effects in the model.

\code{MAsim.real} finally uses existing real or simulated existing data sets to generate simulated data with fixed effect sizes: for each group, the specified number of samples is sampled either with or without replacement from the columns of \code{xdat}; for each gene, the group means are subtracted from the resampled data, so that the groupwise and overall mean for each gene is zero. Then, noise from an appropriate t-distribution is added to each group to break the sum-to-zero constraint in a consistent manner. The specified effect (evenly split between up- and down-regulation) for the differentially expressed genes is again added to the second group.
}
\value{
The functions all return a matrix with \code{ng} rows and \code{n1+n2} columns, except for \code{MAsim.real}, where the default is to return a matrix of the same dimensions as \code{xdat}. The group membership of each column is given by its column name. The matrix has additionally the attribute \code{DE}, which is a logical vector specifying for each gene whether or not it was assumed to be differentially expressed in the simulation.
}
\references{
Smyth G (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3, No. 1, Article 3}
\seealso{\code{\link{EOC}}}
\examples{
# Small examples only
sim1 = MAsim(ng=1000, n=10, p0=0.8)
sim2 = MAsim.var(ng=1000, n1=15, n2=5, p0=0.8)
sim3 = MAsim.smyth(ng=1000, n=10, p0=0.8)

# Assess FDR
eoc1 = EOC(sim1, colnames(sim1), plot=FALSE)
eoc2 = EOC(sim2, colnames(sim2), plot=FALSE)
eoc3 = EOC(sim3, colnames(sim3), plot=FALSE)

# Show
par(mfrow=c(2,2))
plot(eoc1)
plot(eoc2)
plot(eoc3)
OCshow(eoc1, eoc2, eoc3)

# The truth will make you fret
table(eoc1$FDR<0.1, attr(sim1, "DE"))
}
\keyword{models}% at least one, from doc/KEYWORDS

