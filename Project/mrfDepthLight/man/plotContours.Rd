\name{plotContours}
\alias{plotContours}
\title{plotContours}
\description{
PLOTCONTOURS plots the bivariate contours calculated via depthContour.
}
\usage{
plotContours(x, depthContour,data=TRUE)
}
\arguments{
\item{x}{ The (n x 2) data matrix.}
\item{depthContour}{ The result of a call to depthContour. }
\item{data}{Logical. If set to TRUE, the data x is plotted. Defaults to TRUE.}
}
\details{
none
}
\references{
none
}
\seealso{
\code{\link{hdepth}}, \code{\link{depthContour}}
}
\author{P. Segaert.}

\examples{
data(cardata90)
Result = depthContour(x=cardata90, alpha=c(0.02,0.125,0.3) , type="hdepth")
plotContours( x = cardata90, depthContour = Result)
Result = depthContour(x=cardata90, alpha=c(0.1, 0.2, 0.3, 0.4) , type="projdepth")
plotContours( x = cardata90, depthContour = Result)
Result = depthContour(x=cardata90, alpha=c(0.1, 0.2, 0.3, 0.4) , type="sprojdepth")
plotContours( x = cardata90, depthContour = Result)

}

\keyword{Graphical}

