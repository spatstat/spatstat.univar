\name{knots.ewcdf}
\alias{knots.ewcdf}
\title{
  Jump Points of an Empirical Weighted Cumulative Distribution Function
}
\description{
  Extract the knots (jump points) of an empirical weighted cumulative
  distribution function.
}
\usage{
 \method{knots}{ewcdf}(Fn, \dots)
}
\arguments{
  \item{Fn}{
    An empirical weighted cumulative
    distribution function (object of class \code{"ewcdf"}).
}
  \item{\dots}{
    Ignored.
  }
}
\details{
  The function \code{\link[stats]{knots}} is generic. This
  function \code{knots.ewcdf} is the method for
  the class \code{"ewcdf"} of empirical weighted cumulative
  distribution functions. Objects of class \code{"ewcdf"}
  are created by \code{\link{ewcdf}}.

  The jump points (locations of increments) of the function will be
  returned as a numeric vector.
}
\value{
  Numeric vector.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{ewcdf}}, \code{\link{quantile.ewcdf}}
}
\examples{
   x <- c(1, 2, 5)
   w <- runif(3)
   e <- ewcdf(x,w)
   knots(e)
}
\keyword{nonparametric}
\keyword{univar}
