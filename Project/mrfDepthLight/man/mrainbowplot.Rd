\name{mrainbowplot}
\alias{mrainbowplot}
\title{mrainbowplot}
\description{
mrainbowplot plots multivariate data and colors it according to its depth value.
}
\usage{
mrainbowplot(x, depths, col=NULL, LegendTitle="")
}
\arguments{
  \item{x}{An \code{(nx2)} data matrix. \cr
          Observations in the rows and variables in the columns.  }
  \item{depths}{An \code{(nx1)} column vector. \cr
                The depths of the observations in x. The coloring
                is based on these depths.}
  \item{col}{An \code{(mx3)} matrix. \cr
             Colors in rgb format. The coloring of the points is
             based on these values.
             }
  \item{LegendTitle}{Title of the legend.
                    Defaults to \code{"depth"}.
                    }
}
\details{
The plot is made using ggplot2. The plot itself is returned by the function and is fully customisable using standard ggplot2 commands.

}

\author{P. Segaert}
\examples{
data(cardata90)
Result = projdepth(x=cardata90)
plot = mrainbowplot(cardata90, depths = Result$depthZ, LegendTitle = "PD")
library("ggplot2")
plot + ggtitle("Rainbowplot of the cardata using projection depth.")
}

\keyword{Graphical}

