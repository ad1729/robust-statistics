\name{depthContour}
\alias{depthContour}
\title{Depth contours of datasets in dimension \eqn{p \ge 2}.}
\description{
 depthcontour computes depth regions of datasets in dimensions
 \eqn{p \ge 2} for halfspace depth, projection depth or skewness-adjusted
 projection depth. Calculation is exact in two dimensions for the halfspace
 depth. Otherwise the user may specify directions on which to search
 intersections with the depth contour.
}
\usage{
depthContour(x,alpha=NULL, type="hdepth", directions=NULL, options=NULL)
}
\arguments{
	\item{x}{An \eqn{n} by \eqn{p} data matrix.}
	\item{alpha}{A vector containing the depths of which the depth
	             regions have to be computed.}
  \item{type}{The depth used in the computations.
              One of the following options: \code{"hdepth"},
              \code{"projdepth"}, \code{"spdepth"}. \cr
              Defaults to \code{"hdepth"}.}
  \item{directions}{An \eqn{m} by \eqn{p} matrix specifying the directions
                    in which to search intersections. }
  \item{options}{A list of options to pass to hdepth, projdepth or sprojdepth
                 during the approximate algorithm.
                 In addition the following option may be specified:
                  \itemize{
                      \item{ \code{maxIter} \cr
                            The maximum number of iterations in the bisection algorithm used
                            to calculate the depth contour corresponding to level alpha. \cr
                            Defaults to \eqn{100}.
                           }
                          }
                }
}
\details{
 It is first checked whether the data is found to lie in a subspace of
 dimension lower than \eqn{p}. If so, the routine will give a warning, giving
 back the dimension of the subspace together with a direction describing a
 hyperplane containing this subspace.

 In two dimensions for the halfspace depth, dithering is performed when
 the data is not general positionby means of
 random standard Gaussian noise. If after dithering the data is still not
 in general position dithering is repeated on the original undithered
 data.

 In all other cases a biscetion approximation algorithm is performed. On each
 direction a point is searched having depth alpha. Starting limits are obtained
 by projection the data on the direction and taking data point having univariate
 depth alpha. By definition the multivariate depth of this point has to be lower
 or equal to alpha. A second limit is obtained by considering the corresponding
 multivariate location estimate.The maximum number of iterations bisecting the
 current search interval can be specified trought the option argument.
}
\value{
The output consists of list with the following elements for each calculated countour:
	\item{vertices}{The coordinates of the vertices of the depth contour.}
	\item{depth}{The depth of the corresponding depth contour.}
  \item{count}{The number of vertices in the corresponding depth contour.}
	\item{empty}{Logical signaling whether the corresponding depth region is empty. \cr
                \code{FALSE} indicates the depth region is non-empty. \cr
                \code{TRUE} indicates the depth region is empty. \cr
                Only for halfspace depth.}
	\item{dithered}{Logical signaling whether whether dithering has been applied. \cr
                \code{FALSE} indicates no dithering has been performed. \cr
                \code{TRUE} indicates dithering has been applied. \cr
                Only for halfspace depth.}
}
\references{
Ruts, I. and Rousseeuw, P.J. (1996), Computing depth contours of bivariate point clouds. \emph{Computational Statistics and Data Analysis}, \bold{23}, 153-168.

}

\author{P. Segaert based on Fortran code by  P. Rousseeuw , J. Ruts and A. Struyf}

\examples{
data(cardata90)
Result = depthContour(x=cardata90, alpha=0.125, type="hdepth")
plotContours(x=cardata90, depthContour = Result)
}
\keyword{multivariate}
