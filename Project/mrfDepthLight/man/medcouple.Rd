\name{medcouple}
\alias{medcouple}
\title{Medcouple, a robust measure of skewness}
\description{
Computes the medcouple, a robust measure of skewness for univariate data, on each column of a matrix.
}
\usage{
medcouple(x,do.reflect=NULL)
}
\arguments{
  \item{x}{An \eqn{n} by \eqn{p} data matrix.}
  \item{do.reflect}{Logical indicating whether the medcouple should also be              computed on the reflected sample \code{-x}, with final result
  \eqn{(mc}(\code{x})\eqn{-mc(-}\code{x}\eqn{))/2}. \cr
  Defaults \code{TRUE} when \eqn{n<=100} and to \code{FALSE} otherwise.}
}
\details{
The medcouple is a robust measure of skewness yielding values between \eqn{-1} and \eqn{1}. For left- and right-skewed data the medcouple is negative and positive respectively. \cr
The medcouple is defined as the median of the kernel function
\eqn{h(x_i,x_j) = \frac{(x-j - med(x)) - (med(x)-x_i)}{x_j-x_i}}
evaluated over all couples \eqn{(x_i,x_j)} where
\eqn{x_i} is smaller than the median of \code{x}, and \eqn{x_j} is larger than the median of \code{x}. For observations tied to the median, computations are described in the paper by Brys et al. (2004).

To compute the medcouple an algorithm with time complexity \eqn{O(n log n)} is applied. For details, see \url{https://en.wikipedia.org/wiki/Medcouple}.


Missing values are not allowed.
}
\value{
  \item{mc}{A \eqn{p}-vector containing the medcouple of each column of the data matrix \code{x}.}
}
\references{
  Brys G., Hubert M., Struyf A. (2004).
  A robust measure of skewness.
  \emph{Journal of Computational and Graphical Statistics}, \bold{13}, 996--1017.

}
\author{P. Segaert with original code from M. Maechler and G. Brys.}

\examples{
data(bloodfat)
medcouple(bloodfat)
medcouple(bloodfat,do.reflect=TRUE)
}

\keyword{multivariate}
