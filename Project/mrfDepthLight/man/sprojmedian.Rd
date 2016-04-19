\name{sprojmedian}
\alias{sprojmedian}
\title{Location estimates based on projection depth.}
\description{
sprojmedian computes a location estimate of a \eqn{p}-dimensional dataset \code{x} based on the skewness-adjusted projection depths \code{x}.
}
\usage{
sprojmedian(x, sprojDepths=NULL, options=NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{sprojDepths}{Vector containing the skewness-adjusted projection
                  depths of the points in \code{x}.}
  \item{options}{ A list of options to pass to the spdepth routine.
                  See spdepth for more details. }

}
\details{
The algorithm depends on the function \code{spdepth} to calculate the
skew-adjusted projection depths of the dataset \code{x}. If the skewness-adjusted
projection depths of the points have already been calculated they can be
passed as an optional argument to save computation time. If not,
skewness-adjusted projections depths will be calculated and the user can
pass a list with options to the sprojdepth function.
}
\value{
The returned value is a list with the following arguments :
\item{max}{The point of \code{x} with maximal skewness-adjusted projection depth.
           If multiple points have maximum depth, their center
           of gravity is returned}
\item{gravity}{The center of gravity for the 50 percent points with
               maximal skewness-adjusted projection depth. Alternatively called
               the center of gravity of the bag.}
}
\references{
 Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.
 Hubert M, Rousseeuw P., Segaert P. (2015). Multivariate Functional Outlier Detection. \emph{Statistical Methods & Applications}, \bold{X}, X--X.
}
\author{P. Segaert}
\examples{
n<-100
set.seed(123)
x<-matrix(rnorm(n*2),nc=2)
c0<-sprojmedian(x)
}
\keyword{multivariate}
