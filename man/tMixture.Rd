\name{tMixture}
\alias{tMixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a mixture of t-distributions}
\description{
For a vector of individual genewise t-statistics, this functions fits a distribution of central and non-central t-distributions, with the primary goal of estimating the proportion \code{p0} of non-differentially expressed genes.
}
\usage{
tMixture(tstat, n1 = 10, n2 = n1, nq, p0, p1, D, delta, paired = FALSE, 
         tbreak, ext = TRUE, threshold.delta=0.75, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{the vector of genewise t-statistics}
  \item{n1}{number of samples in the first group}
  \item{n2}{number of samples in the second group}
  \item{nq}{the number of components in the mixture that is fitted}
  \item{p0}{a starting value for the proportion of non-differentially expressed genes.}
  \item{p1}{a vector with starting values for the proportions of genes that are differentially expressed with effect size \code{D}.}
  \item{D}{a vector of starting values for the effect sizes of the differentially expressed genes, corresponding to the proportions \code{p1}.}
  \item{delta}{a vector of starting values for the effect sizes of the differentially expressed genes, expressed as non-centrality parameters; this is just a different way of specifying \code{D}, though if both are given, \code{delta} will get priority.}
  \item{paired}{a logical value indicating whether the t-statistics are two-sample or paired.}
  \item{tbreak}{either the number of equally spaced bins for tabulating \code{tstat}, or the explicit break points for the bins, very much like the argument \code{breaks} to function \code{cut}; the default value is the square root of the number of genes.}
  \item{ext}{a logical value indicating whether to extend the bins, i.e. to set the lowest bin limit to -infinity and the largest bin limit to inifinity.}
  \item{threshold.delta}{mixture components with an estimated absolute non-centrality parameter \code{delta} below this value are considered to be too small for independent estimation; these components and their corresponding \code{p1} are pooled with the null-component and \code{p0}, see Details.}
  \item{\dots}{additional arguments that are passed to \code{optim} to control the optimization.}
}
\details{
The minimum parameter that needs to be specified is \code{nq} - if nothing else is given, the proportions are equally distributed between \code{p0} and the \code{p1}, and the noncentrality parameters are set up symmetrically around zero, e.g. \code{nq=5} leads to equal proportions of 0.2 and noncentrality parameters -2, -1, 1, and 2. If any of \code{p1}, \code{D}, or \code{delta} is specified, \code{nq} is redundant and will be ignored (with a warning). \code{tMixture} will in general make a valiant effort to deduce valid starting values from any combination of \code{nq}, \code{p0}, \code{p1}, \code{D}, and \code{delta} specified by the user, and will complain if that is not possible. 

The fitting problem that this function tries to solve is badly conditioned, and will in general depend on the precise set of starting values. Multiple runs from different starting values are usually a good idea. We have found however, that the model seems fairly robust towards misspecification of the number of components, at least when estimating \code{p0}. What happens when too many components are specified is that some of the nominally noncentral t-distributions describing the behaviour of differentially expressed genes are fitted with noncentrality parameters very close to zero, and the true \code{p0} gets spread out between the nominal \code{p0} and the almost-central components. Adding up these different contributions usually gives a similar solution to re-fitting the model with fewer components. The cutoff for the size of  non-centrality parameters that can be estimated realistically is specified via \code{threshold.delta}, whose default value is based on a small simulation study reported in Pawitan et al. (2005); see Examples. (Note that the AIC can also be helpful in determining the number of components.)  
}
\value{
A list with the following components:
  \item{p0.est}{the estimated proportion of non-differentially expressed genes, after collapsing components with estimated non-centrality sizes below \code{threshold.delta}.}
  \item{p0.raw}{the estimated proportion before collapsing the components.} 
  \item{p1}{the estimated proportions of differentially expressed genes corresponding to the effect sizes, relating to \code{p0.raw}.}
  \item{D}{effect sizes of the differentially expressed genes in multiples of the gene-by-gene standard deviation.}
  \item{delta}{effect sizes of the differentially expressed genes expressed as the noncentrality parameter of the corresponding noncentral t-distribution.}
  \item{AIC}{the AIC value for the maximum likelihood fit.}
  \item{opt}{The output from \code{optim}, giving details about the optimization process.}
}
\references{
Pawitan Y, Krishna Murthy KR, Michiels S, Ploner A (2005) Bias in the estimation of false discovery rate in microarray studies, Bioinformatics.
}
\author{Y. Pawitan and A. Ploner}
\seealso{\code{\link{tstatistics}}, \code{\link{EOC}}, \code{\link{optim}}}
\examples{
# We simulate a small example with 5 percent regulated genes and
# a rather large effect size
set.seed(2011)
xdat = matrix(rnorm(50000), nrow=1000)
xdat[1:25, 1:25] = xdat[1:25, 1:25] - 2
xdat[26:50, 1:25] = xdat[26:50, 1:25] + 2
grp = rep(c("Sample A","Sample B"), c(25,25))
# Use a helper function for the test statistics
myt = tstatistics(xdat, grp)$tstat
r1 = tMixture(myt, n1=25, nq=3)
r1

# Equivalently, we could have specified the same set of starting values 
# as follows:
# r1 = tMixture(myt, n1=25, p0=1/3, p1=c(1/3, 1/3), delta=c(-1,1))

# Alternative starting value for p0, other starting values are filled in
r2 = tMixture(myt, n1=25, nq=3, p0=0.80)
r2

# Specification of too many components usually leads to spurious
# noncentral components like here - note the difference between
# p0.est and p0.raw!
r3 = tMixture(myt, n1=25, nq=5)
r3            

# We simulate a data in a paired setting
# Note the arrangement of the columns
set.seed(2012)
xdat = matrix(rnorm(50000), nrow=1000)
ndx1 = seq(1,50, by=2)
xdat[1:25, ndx1] = xdat[1:25, ndx1] - 2
xdat[26:50, ndx1] = xdat[26:50, ndx1] + 2
grp = rep(c("Sample A","Sample B"), 25)
# Use a helper function for the test statistics
myt = tstatistics(xdat, grp, paired=TRUE)$tstat
r1p = tMixture(myt, n1=25, nq=3)
r1p
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
