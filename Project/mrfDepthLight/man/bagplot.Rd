\name{bagplot}
\alias{bagplot}
\title{Draws a bagplot of bivariate data.}
\description{
This function draws a bagplot of bivariate data, based on the result of a call to \code{compBagplot}. The bagplot is a generalisation of the univariate boxplot to bivariate data.
}
\usage{
bagplot(CompBagResult,colorbag=NULL,colorloop=NULL,colorchull=NULL,
        databag=TRUE,dataloop=TRUE,plotFence=FALSE)
}
\arguments{
  \item{CompBagResult}{The return of a call to \code{compBagplot}.}
  \item{colorbag}{The color of the bag (which contains the 50\% observations with largest depth).}
	\item{colorloop}{The color of the loop (which contains the regular observations).}
  \item{colorchull}{When the bagplot is based on halfspace depth the depth region with maximal depth is plotted. This argument controls its color.}
  \item{databag}{Logical indicating whether data points inside the bag need to be plotted. \cr
  Defaults to \code{TRUE}.}
  \item{dataloop}{Logical indicating whether data points inside the fence need to be plotted. \cr
  Defaults to \code{TRUE}.}
    \item{plotFence}{Logical indicating whether the fence should be plotted. \cr
  Defaults to \code{FALSE}.}
}
\details{
The bagplot has been proposed by Rousseeuw et al. (1999) as a generalisation of the boxplot to bivariate data. It is constructed based on halfspace depth. Similar graphical representations can be obtained by means of other depth functions, as illustrated in Hubert and Van der Veeken (2008). See \code{\link{compBagplot}} for details.

The deepest point is indicated with a "+" sign, the outlying observations with red stars.

The plot is made using \code{ggplot2}. The plot itself is returned by the function and is fully customisable using standard \code{ggplot2} commands.
}
\seealso{
  \code{\link{compBagplot}}, \code{\link{hdepth}}, \code{\link{projdepth}}, \code{\link{sprojdepth}}.
}
\references{
Rousseeuw, P.J., Ruts, I., Tukey, J.W. (1999). The bagplot: a bivariate boxplot. \emph{The American Statistician},\bold{53}, 382--387.

Hubert, M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.

Hubert, M. Rejoinder ??

}

\author{P. Segaert}

\examples{
data(bloodfat)
# Bagplot based on halfspace depth
Result = compBagplot(bloodfat)
bagplot(Result)
bagplot(Result,databag=FALSE,dataloop=FALSE)
bagplot(Result,colorbag=rgb(0.2,0.2,0.2),colorloop="green")
# Bagplot based on projection depth
Result = compBagplot(bloodfat, type="projdepth")
bagplot(Result)

# The output of blagplot is a ggplot2 object which can be modified.
library(ggplot2)
Plot = bagplot(Result)
Plot + xlab("Your personalized x- label") + ggtitle("Your own title. ")

}
\keyword{Graphical}
