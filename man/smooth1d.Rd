\name{smooth1d}
\alias{smooth1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothing a vector of counts}
\description{
This function takes a vector of counts and uses a mixed model approach to smooth it. A common use of this is smoothing binned counts of an observed quantity prior to estimating its density nonparametrically through the relative frequencies. 
}
\usage{
smooth1d(y, sv2 = 0.1, err = 0.01, verb = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the vector of counts}
  \item{sv2}{the user-specified starting value for the variance of the random effects, see Details.}
  \item{err}{Tolerance for convergence, see Details}
  \item{verb}{logical value indicating whether to print diagnostics.}
}
\details{
The smoothing assumes that the counts are Poisson from a generalized linear mixed model, where the second differences are normally distributed. Using the extended likelihood approach described in Pawitan (2001) and the initial estimate \code{sv2} for the variance of the random effects, the routine iteratetively optimizes the fixed and random contributions to the extended likelihood, until the estimate for the variance convergences with tolerance \code{err}. The result is quite stable within a reasonable range of starting values and tolerances, and the function can be used for fairly automatic smoothing ((i.e. withou fixing a bandwidth parameter).
}
\value{
A list with three components:
\item{fit}{the smoothed counts}
\item{df}{the degrees of freedom used for smoothing at convergence}
\item{sv2}{the estimated variance at convergence, equivalent to \code{df}.}
}
\references{
Pawitan Y.(2001) \emph{In All Likelihood}, Oxford University Press, ch. 18.11
}
\author{Y. Pawitan and A. Ploner}
\seealso{\code{\link{fdr1d}}}
\examples{
# Stupid dummies, obviously
smooth1d(1:10)
smooth1d(1:10, sv2=1)
}
\keyword{smooth}% at least one, from doc/KEYWORDS

