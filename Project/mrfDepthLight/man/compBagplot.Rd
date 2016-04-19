\name{compBagplot}
\alias{compBagplot}
\title{Performs all computations to draw a bagplot}
\description{
Computes all elements constituting a bagplot, a generalisation
of the univariate boxplot to bivariate data. The bagplot can be
computed based on halfspace depth, projection depth or skewness-adjusted projection depth.
}
\usage{
compBagplot(x,type="hdepth", sizesubset=500, extraDirections = FALSE, options=NULL)
}
\arguments{
	\item{x}{An \eqn{n} by 2 data matrix with observations in the rows and variables in the columns.}

  \item{type}{Determines the depth function used to draw the bagplot: \code{"hdepth"} for halfspace depth, \code{"projdepth"} for projection depth and \code{"spdepth"} for skewness-adjusted projection depth}. \cr
              Defaults to \code{"hdepth"}.

	\item{sizesubset}{When computing the bagplot based on the halfspace depth,
                    the size of the subset used to perform the main
                    computations. \cr
                    Defaults to \eqn{500}.
                    }
  \item{extraDirections}{Logical indicating whether additional directions should
                         be considered in the calculation of the fence for the
                         bagplot based on projection depth or skewness-adjusted
                         projection depth.\cr
                         Defaults to \code{FALSE}.
                        }
  \item{options}{A list of options to pass to the projdepth or spdepth routine.
                 In addition the following option may be specified:
                  \itemize{
                      \item{ \code{maxIter} \cr
                            The maximum number of iterations in the bisection algorithm used
                            to calculate the depth contour corresponding to the cutoff \cr
                            Defaults to \eqn{100}.
                           }
                          }

                }

}
\details{
The bagplot has been proposed by Rousseeuw et al. (1999) as a generalisation of the boxplot to bivariate data.  It is constructed based on halfspace depth. In the original format the deepest point is indicated by a "+" and it is contained in the bag which is the depth region containing the 50\% observations with largest depth. The fence is obtained by inflating the bag (relative to the deepest point) by a factor three. The loop is the convex hull of the observations of \code{x} inside the fence. Observations outside the fence are flagged as outliers and plotted with a red star. This function only computes all the components constituting the bagplot. The bagplot itself can be drawn using the \code{bagplot} function.

However the bagplot may also be defined using other depth functions. When using projection depth or skewness-adjusted projection depth the bagplot is build as follows. The center corresponds to the observation with depth. The bag is constructed as the convex hull of the fifty percent points with largest depth. Outliers are identified as points with depth smaller than a cutoff value.
This cutoff value corresponds to the depth of point with outlyingness value equalling the square root of \eqn{.99} quantile of a chi-square distirbution with 2 degrees of freedom multiplied by the median of the outlyingness values of the depths of the data matrix \code{x}. The fence is calculated as the convex hull of the non-outlying points. The loop is approximated by searching on the rays from the center trough the cornerpoints of the bag and the fence for the point having depth equal to the cutoff depth. For a better approximation the used set the flag \code{ExtraDirections} to \code{TRUE} such that an additional 250 equaly spaced directions on the circle are considered.

The calculation of the bagplot based on halfspace depth can be time
 consuming. Therefore it is possible to limit the bulk of the calculations
 on a random subset of the data. Calculation of the halfspace median and
 the bag is based on this random subset. The number of points in this
 subset can be controlled by the optional argument \code{sizesubset}.
If projection depth or the skewness-adjusted projection depth is used to construct the bagplot, the outliers are defined as those observations whose outlyingness or adjusted outlyingness exceed the corresponding cutoff value.

 It is first checked whether the data is found to lie in a subspace of
 dimension lower than two. If so, the routine will give a warning, giving
 back the dimension of the subspace together with a direction describing a
 hyperplane containing this subspace.



}
\value{
 A  list with components: \cr
	\item{center}{Center of the data. \cr
                       When \code{type=halfspace}, this corresponds with the
                       Tukey median.
                       In other cases this point corresponds to the point
                       with maximal depth.
	}
  \item{chull}{When \code{type=halfspace}, these are the vertices of
               the region with maximal halfspace depth.
               In other cases this is a null vector.
  }
  \item{bag}{The points on the contour of the bag.}
  \item{datatype}{An \code{(nx3)} matrix with the third column indicating the
                 position of the point in the bagplot. 2 for observations in
                 the bag, 1 for the observations in the fence and three for
                 outliers. \cr
                 Note that points may not be in the same order as in \code{x}.

  }
	\item{flag}{0 for outliers and 1 for regular points.}
  \item{depth}{The depth of the observations in the dataset.}
}
\seealso{
  \code{\link{bagplot}}, \code{\link{hdepth}}, \code{\link{projdepth}}, \code{\link{sprojdepth}}.
}
\references{
Rousseeuw, P., Ruts, I., and Tukey, J.W. (1999). The bagplot: A bivariate boxplot. \emph{The American Statistician},\bold{53}, 382--387.

Hubert, M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.
}

\author{P. Segaert based on Fortran code by P. Rousseeuw, I. Ruts and A. Struyf.}

\examples{
data(bloodfat)
Result = compBagplot(bloodfat);
bagplot(Result)
}
\keyword{multivariate}
