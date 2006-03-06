\name{summary.fdr.result}
\alias{summary.fdr.result}
\alias{p0}
\alias{topDE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Display functions for local fdr output}
\description{
Display functions for output from \code{fdr1d} and \code{fdr2d}, summarizing the output, displaying the proportion of non-differentially expressed genes and extracting the list of top-regulated genes. 
}
\usage{
summary.fdr.result(object, ...)

p0(x, how = FALSE)

topDE(x, co = 0.1)
}
\arguments{
  \item{x, object}{an object of class \code{fdr.result} created by \code{fdr1d} or \code{fdr2d}.}
  \item{how}{a logical value indicating whether to return only the numerical value of the proportion of non-differentially expressed genes, or a list whose second element indicates whether the proportion was estimated from the data or supplied by the user.}
  \item{...}{extra arguments, currently unused}
  \item{co}{cutoff for either FDR or fdr}
}
\value{
For \code{summary.fdr.result}, a list with the summary items.

For \code{p0}, either a numerical value or a list with two elements, depending on the value of parameter \code{how}.

For \code{topDE}, the genes that have FDR (\code{EOC}) or fdr (\code{fdr1d}, \code{fdr2d}) less or equal than \code{co}, sorted by FDR or fdr respectively.
}
\author{A. Ploner}
\seealso{\code{\link{fdr1d}}, \code{\link{fdr2d}}, \code{\link{EOC}}}
\examples{
# Create object res1d
example(fdr1d)

summary(res1d)
p0(res1d)
p0(res1d, how=TRUE)
topDE(res1d)
}
\keyword{print}% at least one, from doc/KEYWORDS
\keyword{utilities}% at least one, from doc/KEYWORDS

