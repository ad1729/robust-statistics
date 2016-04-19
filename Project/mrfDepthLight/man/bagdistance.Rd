\name{bagdistance}
\alias{bagdistance}
\title{Bagdistance}
\description{
Computes the bagdistance of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. 
}
\usage{
bagdistance(x, z = NULL, options = list())
}
\arguments{
	\item{x}{An \eqn{n} by \eqn{p} data matrix with observations
	         in the rows and variables in the columns.}
  \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
           rowwise the points \eqn{z_i} for which to compute
           the bagdistance. If \code{z} is not specified,
           it is set equal to \code{x}.}
  \item{options}{A list of available options:
                \itemize{
                      \item{ \code{maxIter} \cr
                            The maximum number of steps in the bisection algorithm to find the intersection point \eqn{c_x} (see Details).\cr
                            Defaults to \eqn{100}.
                            }
                \item{}
                All options may be specified that are passed to the \code{hdepth} function, see \code{hdepth} for details. Note that the option parameter \code{approx} is by default set to \code{TRUE} to save computation time.}}
}
\value{
  A list with components: \cr
  \item{bagdistance}{The bagdistance of the points of \code{z} with respect to the data matrix \code{x}.}
  \item{converged}{Vector of length \code{m} indicating for each point of \code{z} whether the bisection algorithm converged within the maximum number of steps specified by \code{maxIter} in the \code{options} list.}
  \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}
  \item{type}{bagdistance}
}
\details{
The bagdistance has been introduced in Hubert et al. (2015). 
To compute the bagdistance of a point \eqn{z_i}, we first consider the bag of \code{x} which is the 50\% central halfspace depth region. This is the smallest depth region that contains at least half of the observations. Next, we consider the ray from the halfspace median \eqn{\theta} through \eqn{z_i}, and we define \eqn{c_x} as the intersection of this ray and the boundary of the bag. The bagdistance of \eqn{z_i} to \code{x} is then given by the ratio between
the Euclidean distance of \eqn{z_i} to the halfspace median and the Euclidean distance of \eqn{c_x} to the halfspace median.
The bagdistance does not assume symmetry and is affine invariant.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it. Otherwise, the function first computes the halfspace depth and the halfspace median of \code{x}. Additional options may be passed to the \code{hdepth} routine by specifying them in the \code{option} list argument.

Depending on the dimensions different algorithms are used. For \eqn{p=1} and \eqn{p=2} the bagdistance is computed exactly, based on the exact algoritm to compute the contours of the bag (see the \code{depthContour} function). In dimensie enkel exact als options$approx=FALSE, zoals gevraagd in de TODO lijst.

For \eqn{p \ge 3}, the intersection point \eqn{c_x} is approximated by searching on each ray the point whose depth is equal to the median of the depth values of \code{x}. As the halfspace depth is monotone decreasing along the ray, a bisection algorithm is used. To save computation time, this algorithm may also be applied in the setting \eqn{p=2} by setting \code{options$approx=TRUE}. 

To flag outliers,  ...

}\seealso{
  \code{\link{depthContour}}, \code{\link{hdepth}}, \code{\link{bagplot}}
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection. \emph{Statistical Methods & Applications}, \bold{24}, 177--202.

}

\author{P. Segaert. }

\examples{
# Bivariate dataset
data(cardata90)
# For the matrix \code{z}, we generate 100 random points and shift them to the
# halfspace median of the car dataset. 
TestData <- matrix(rnorm(200,sd=50),ncol=2)
TestData <- sweep(TestData, MARGIN=2, FUN="+", hdepthmedian(cardata90)$median)
Distances <- bagdistance(x = cardata90,z = TestData)
}

\keyword{multivariate}

