\name{hdepthmedian}
\alias{hdepthmedian}
\title{Location of maximal halfspace depth}
\description{hdepthmedian computes the halfspace median and its corresponding
 halfspace depth for a \eqn{p}-dimensional dataset \code{x}.  Calculation is exact when \eqn{p \le 2} and approximate when \eqn{p > 2}.
}
\usage{
hdepthmedian(x,maxdir=NULL)
}
\arguments{
  \item{x}{An \eqn{n} by \eqn{p} data matrix .
          Observations in the rows and variables in the columns.}
	\item{maxdir}{The number of projection used in the approximate algorithm. \cr
                Defaults to \eqn{250p}. \cr
                Approximate algorithm only. }
}
\details{
 It is first checked whether the data is found to lie in a subspace of
 dimension lower than \eqn{p}. If so, the routine will give a warning, giving
 back the dimension of the subspace together with a direction describing a
 hyperplane containing this subspace.

 For the case \eqn{p=2} dithering is performed by means of random
 standard Gaussian noise when the data is not in general position.
 If after dithering the data is still not in general position
 dithering is repeated on the original undithered data.
 The approximate algorithm uses an iterative projection procedure.

}
\value{
 A  list containing: \cr
  \item{median}{The coordinates of the halfspace median.\cr
            Approximated when \eqn{p>2}}
  \item{depth}{The depth of the halfspace median. \cr
            Approximated when \eqn{p>2}}
  \item{dithered}{Logical signaling whether dithering has been applied. \cr
                  \code{FALSE} indicates no dithering has been performed. \cr
                  \code{TRUE} indicates dithering has been applied. \cr}
	\item{ndir}{The number of projections used by the algorithm. Due to the possibility of singularity of certain \eqn{p} subsamples it is possible not all maxdir directions are evaluated. \cr
  Approximate algorithm only.}

	\item{AlgStopFlag}{Indicates which stopping rule is used by the algorithm.\cr
                     0 indicates the maximum number of projections was reached \cr
                     1 indicates no improvement of the location estimate was
                     made after ntry steps.\cr
                    Approximate algorithm only.}
}
\references{
Rousseeuw P., Ruts I. (1998), Constructing the bivariate Tukey median, \emph{Statistica Sinica}, \bold{8}, 827--839.

Struyf, A. and Rousseeuw, P. (2000), High-dimensional computation of the deepest location, \emph{Computational Statistics and Data Analysis}, \bold{34}, 415--436.
}

\author{K. Vakili and P. Segaert based on Fortran code by P. Rousseeuw, J. Ruts and A. Struyf}

\examples{
data(cardata90)
Result = hdepthmedian(x=cardata90)
}
\keyword{multivariate}
