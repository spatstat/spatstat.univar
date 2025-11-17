\name{knots.ewcdf}
\alias{knots.ewcdf}
\alias{knots.ecdf}
\title{
  Jump Points of an Empirical Weighted Cumulative Distribution Function
}
\description{
  Extract the knots (jump points) of an empirical cumulative
  distribution function.
}
\usage{
 \method{knots}{ewcdf}(Fn, \dots)
 \method{knots}{ecdf}(Fn, \dots)
}
\arguments{
  \item{Fn}{
    An empirical cumulative distribution function
    (object of class \code{"ecdf"} or \code{"ewcdf"}).
}
  \item{\dots}{
    Ignored.
  }
}
\details{
  The function \code{\link[stats]{knots}} is generic.

  The function \code{knots.ecdf} is the method for the class \code{"ecdf"}
  of empirical cumulative distribution functions; objects of this class
  are created by \code{\link[stats]{ecdf}}).

  The function \code{knots.ewcdf} is the method for
  the class \code{"ewcdf"} of empirical weighted cumulative
  distribution functions. Objects of class \code{"ewcdf"}
  are created by \code{\link{ewcdf}}.

  The jump points (locations of increments) of the function \code{Fn}
  will be returned as a numeric vector.
}
\value{
  Numeric vector.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[stats]{ecdf}}, 
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
