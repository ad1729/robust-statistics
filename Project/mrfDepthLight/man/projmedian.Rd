\name{projmedian}
\alias{projmedian}
\title{Location estimates based on projection depth.}
\description{
projmedian computes a location estimate of a \eqn{p}-dimensional dataset \code{x} based on the projection depths \code{x}.
}
\usage{
projmedian(x,projDepths=NULL, options=NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{projDepths}{Vector containing the projection depths of the points in \code{x}.}
  \item{options}{ A list of options to pass to the projdepth routine.
                  See projdepth for more details. }

}
\details{
The algorithm depends on the function projdepth to calculate the
 projection depths of the dataset \code{x}. If the projection depths of the
 points have already been calculated they can be passed as an optional
 argument to save computation time. If not, projections depths will be
 calculated and the user can pass a list with options to the
 projdepth function.

 It is first checked whether the data is found to lie in a subspace of
 dimension smaller than \code{p}. If so, the routine will give a warning, as well
 as the dimension of the subspace and a direction which is orthogonal to a
 hyperplane containing this subspace.

}
\value{
The returned value is a list with the following arguments :
\item{max}{The point of x with maximal projection depth.
           If multiple points have maximum depth, their center
           of gravity is returned}
\item{gravity}{The center of gravity for the 50 percent points with
               maximal projection depth. Alternatively called the
               center of gravity of the bag.}
\item{Huber}{The weighted center of gravity of all points.
                      The weights are defined by the Huber
                      function with parameter
                          alpha = 1/(1+sqrt(qchisq(0.95,p)))
                      Observations for which the projection depth
                      is more than alpha receive weight 1,
                      other points receive weight
                      (sqrt(qchisq(0.95,p))/outlyingness)^2.}
}
\references{
 Zuo. Y (2003). Projection based depth functions and associated medians.  \emph{The Annals of Statistics}, \bold{31}, 1460--1490.
}
\author{P. Segaert}
\examples{
data(bloodfat)
Result = projmedian(bloodfat)
}
\keyword{multivariate}
