\name{plot.fdr2d.result}
\alias{plot.fdr2d.result}
\alias{Volcanoplot}
\alias{Tornadoplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting the bivariate local false discovery results}
\description{
These functions provide different ways of plotting the output \code{fdr2d}.
}
\usage{
\method{plot}{fdr2d.result}(x, levels, nr.plot = 20, add = FALSE, grid = FALSE, 
  pch = ".",xlab, ylab, vfont = c("sans serif", "plain"), lcol = "black", ...)

Tornadoplot(x, levels, nr.plot = 20, label = FALSE, constrain = FALSE, 
  pch = ".", xlab, ylab, vfont = c("sans serif", "plain"), lcol = "black", ...)
  
Volcanoplot(x, df, levels, nr.plot = 20, label = FALSE, constrain = FALSE, 
  pch = ".", xlab, ylab, vfont = c("sans serif", "plain"), lcol = "black", ...)	   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object created by \code{fdr2d}.}
  \item{df}{the appropriate degrees of freedom for a two-sample t-test with equal variances (in order to provide p-values for the volcano plot).}  
  \item{levels}{vector of levels for drawing fdr isolines}
  \item{nr.plot}{number of equidistant breaks defining a two-dimensional grid for smoothing isolines, see Details.}
  \item{add}{logical value indicating whether to create a new plot, or to add to an existing plot.}
  \item{grid}{logical value indicating whether the original grid used for estimating the local fdr should be shown.}
  \item{label}{logical value indicating whether to draw labels on the isolines.}
  \item{constrain}{logical flag indicating whether transformed fdr values should be made monotnously decreasing with the absolute size of the first test statistic, see \code{fdr2d}.}
  \item{pch, xlab, ylab}{the usual graphical parameters}
  \item{vfont}{vector font specification for labelling the isolines, see \code{contour}.}
  \item{lcol}{colour used for drawing the isolines}
  \item{\dots}{extra graphical parameters passed to \code{plot.default}.}
}
\details{
The plot format is basically a scatter plot of the observed test statistics, overlayed with isolines showing the estimated fdr. The generic plot function displays the original test statistics that are used to estimate the fdr, i.e. the two-sample t-statistics and the logarithmized standard errors; the other plots use different, but mathematically equivalent test statistics:
\itemize{
\item mean difference and logarithmized standard error for the tornado plot,
\item mean difference and \eqn{-\log10(p)} for the volcano plot, where \eqn{p} is the p-value from the standard two-sample t-test.
}

By default, the estimated fdr isolines are smoothed and cropped to the convex hull of the observed test statistics by using \pkg{akima}. This is entirely a graphical pre-processing step which produces smoother isolines and enforces sanity at the edges of the observed distribution; it does not change the estimated fdr at all. This graphical smoothing is controlled via the argument \code{nr.plot}, which specifies the grid size, with lower values resulting in stronger smoothing. In order to suppress graphical smoothing, set \code{nr.plot} to zero.

Note that the test statistics and the fdr for the volcano- and tornado plots are not computed from scratch, but rather through transformation of the original results. Specifically, the isolines in these plots are also transformed; this has the unfortunate side effect that the labelling of isolines in these plots is not nearly as pretty as the standard provided by \code{contour}. This  functionality is currently not available outside of \code{contour}, and our implementation in \code{DrawContourlines} frankly leaves a lot to be desired. We apologize for the inconvenience.}    

\value{
The original \code{x}, invisibly.
}
\references{
Ploner A, Calza S, Gusnanto A, Pawitan Y (2005) Multidimensional local false discovery rate for micorarray studies. \emph{Submitted Manuscript}.
}
\author{A. Ploner}
\seealso{\code{\link{fdr2d}}, \code{\link{DrawContourlines}}}
\examples{
# Create res2d
example(fdr2d)

par(mfrow=c(2,2))
plot(res2d, main="Generic plotting")
Volcanoplot(res2d, df=length(grp)-2, main="Volcano plot", label=TRUE)
Tornadoplot(res2d, main="Tornado plot", label=TRUE)
# This is without graphical smoothing
plot(res2d, main="Generic plotting, raw", nr.plot=0)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{aplot}% __ONLY ONE__ keyword per line
