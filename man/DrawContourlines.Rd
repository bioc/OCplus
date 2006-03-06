\name{DrawContourlines}
\alias{DrawContourlines}
\title{Draw and label a set of pre-calculated isolines}
\description{
This function  draws and labels isolines computed by \code{contourLines}, though the labelling is done very clumsily and with a specialized application in mind.
}
\usage{
DrawContourlines(x, label = FALSE, cex = 0.7, vfont = c("sans serif", "bold"), ...)
}
\arguments{
  \item{x}{a list of isolines as produced by \code{contourLines}.}
  \item{label}{a logical value indicating whether to label the isolines.}
  \item{cex}{size of labels}
  \item{vfont}{a vector font specification for the labels as in \code{contour}.}
  \item{\dots}{extra arguments to \code{lines}}
}
\details{
This routine is used by \code{Tornadoplot} and \code{Volcanoplot} to draw and label isolines that were computed via contourLines and afterwards transformed. The problem is that all the nice options that \code{contour} has for labelling isolines are not avaiable independently, so this function uses the following crude procedure that kind of works for the intended applications: 
\itemize{
\item{isoline completely left of zero}{label the leftmost point;}
\item{isoline completely right of zero}{label the rightmost point;}
\item{isoline crosses zero horizontally}{label the topmost point.}
}
Hopefully, one of these days someone will come up with a nice general-purpose function for doing all the nifty stuff that \code{contour} offers.
}
\author{A. Ploner}
\seealso{\code{\link{contour}}, \code{\link{contourLines}}, \code{\link{Tornadoplot}}}
\keyword{aplot}% at least one, from doc/KEYWORDS

