\name{plot.FDR.result}
\alias{plot.FDR.result}
\title{Plot the empirical FDR as a function of the cutoff level}
\description{
Plots the output from \code{EOC}. The resulting graph is the empirical counterpart to those produced by \code{TOC}, i.e. the estimated FDR as a function of the cutoff-level on the t-statistic.
}
\usage{
plot.FDR.result(x, add=FALSE, sensitivity.show = TRUE, legend.show = FALSE, 
     xlim, ylim = c(0, 1), xlab, ylab, main, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object created by \code{EOC}}
  \item{add}{logical value indicating whether to add to an existing plot or start a new one} 
  \item{sensitivity.show}{logical value indicating whether to show the classical sensitivity for testing one hypothesis as a function of the cutoff level.}
  \item{legend.show}{logical value indicating whether to add a legend to the plot}
  \item{xlim, ylim}{limits for the horizontal and vertical axis}
  \item{xlab, ylab}{axis labels}
  \item{main}{plot title}
  \item{\dots}{the usual graphical parameters, passed to \code{plot}}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
%\references{ ~put references to the literature/web site here ~ }
\author{A. Ploner}
%\note{ ~~further notes~~ }
\seealso{\code{\link{EOC}}}
\examples{
# We simulate a small example with 5 percent regulated genes and
# a rather large effect size
set.seed(2003)
xdat = matrix(rnorm(50000), nrow=1000)
xdat[1:25, 1:25] = xdat[1:25, 1:25] - 2
xdat[26:50, 1:25] = xdat[26:50, 1:25] + 2
grp = rep(c("Sample A","Sample B"), c(25,25))

# Compute the EOC without plotting
ret = EOC(xdat, grp, plot=FALSE) 

# Some possible arrangements
par(mfrow=c(2,2))
plot(ret)
plot(ret, legend=TRUE)
plot(ret, sensitivity=FALSE)

}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{aplot}% __ONLY ONE__ keyword per line
