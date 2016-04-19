\name{cardata90}
\alias{cardata90}
\title{Cardata from Consumer Reports in 1990.}
\description{
  Subset form data on cars taken from pages 235-255, 281-285 and 287-288 of the april 1990 Consumer Reports Magazine.
}
\usage{
  data(cardata90)
}
\format{
  Following variables are included:
    \describe{
      \item{\code{Weight}}{ Weight of the car }
      \item{\code{Disp}}{ Engine displacement in cubic inches. }
    }
}
\source{
\emph{Consumer Reports}, April 1990, pp. 235-288.

Chambers, J.M. and Hastie T.J. (1993). Statistical models in S \emph{Londen: Chapman and Hall}, 46-47.
}
\examples{
  data(cardata90)
  plot(cardata90)
}
\keyword{datasets}
