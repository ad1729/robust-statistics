\name{sprojdepth}
\alias{sprojdepth}
\title{Skewness-adjusted projection depth}
\description{
 sprojdepth computes the skewness-adjusted projection depth of \eqn{p}-dimensional
 points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}.
 Skewness-adjusted projection depth is based on the adjusted
 outlyingness (see adjOutlyingness) and is computed as \eqn{1/(1+AO)}.
}
\usage{
sprojdepth(x,z=NULL,options=NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
    \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
             rowwise the points \eqn{z_i} for which to compute
             the projection depth. If \code{z} is not specified,
             it is set equal to \code{x}.}
  \item{options}{ A list of options to pass to the adjOutlyingness routine.
                  See adjOutlyingness for the full list of options. }
}
\details{
 It is first checked whether the data is found to lie in a subspace of
 dimension lower than \eqn{p}. If so, the routine will give a warning, giving
 back the dimension of the subspace together with a direction describing a
 hyperplane containing this subspace.
}
\value{
  A list with components: \cr
      \item{depthX}{Vector of length \eqn{n} giving the
                    skewness-adjusted projection depth of
                    the observations in x.}
      \item{depthZ}{Vector of length \eqn{m} giving the
                    skewness-adjusted projection depth of
                    the points in z.}
      \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"},
                             the number of \eqn{p}-subsets that span a subspace of
                             dimension smaller than \eqn{p-1}. In that case the orthogonal
                             direction can not be uniquely determined. This is an indication
                             that the data are not in general position.
                             When the input parameter type is equal to \code{"Rotation"} it
                             is possible that two randomly selected points of the data coincide
                             due to ties in the data. In this case this value signals how many
                             times this is the case.
                             }
  \item{cutoff}{Points whose skew-adjusted projection depth is smaller
                than this cutoff can be considered as outliers.}
  \item{flagX}{Observations of \code{x} whose adjusted outlyingness exceeds the cutoff receive a flag \code{FALSE}, regular observations receive a flag \code{TRUE}.}
  \item{flagZ}{Points of \code{z} whose adjusted outlyingness exceeds the cutoff receive a flag equal to \code{FALSE}, otherwise they receive a flag \code{TRUE}.}
  \item{dimension}{When the data \code{x} are lying in a
                   lower dimensional subspace, the dimension
                   of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a
                    lower dimensional subspace, a direction orthogonal
                    to this subspace. When a direction is found with
                    \code{IQR = 0} this equals to the direction in question.}
  \item{SubspaceIndicator}{When a direction is found for which the IQR
                           after projection is zero this value is a vector
                           of length \eqn{n} signaling which of the points
                           are different from the median after projection
                           on a certain directions. These points receive an
                           indicator 1. Only returned when such a direction
                           is found.}
}
\seealso{
  \code{\link{adjOutlyingness}}.
}
\references{
 Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.
 Hubert M, Rousseeuw P., Segaert P. (2015). Multivariate Functional Outlier Detection. \emph{Statistical Methods & Applications}, \bold{X}, X--X.
}
\author{P. Segaert with original code from M. Maechler, G. Brys, K. Vakili}
\examples{
data(bloodfat)
Result = sprojdepth(bloodfat)
mrainbowplot(bloodfat, depth = Result$depthX)

}
\keyword{multivariate}
