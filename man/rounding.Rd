\name{rounding} 
\alias{rounding}
\alias{rounding.default}
\title{
  Detect Numerical Rounding 
}
\description{
  Given a numeric vector,
  determine whether the values have been rounded
  to a certain number of decimal places.
}
\usage{
rounding(x)

\method{rounding}{default}(x)
}
\arguments{
  \item{x}{
    A numeric vector, or an object containing numeric spatial coordinates.
  }
}
\details{
  The function \code{rounding} is generic.
  Its purpose is to determine whether numerical values have been rounded
  to a certain number of decimal places.

  The \pkg{spatstat} family of packages provides methods for
  \code{rounding} for various spatial objects.
  
  For a numeric vector \code{x}, the default method \code{rounding.default}
  determines whether the values in \code{x} have been rounded
  to a certain number of decimal places.
  \itemize{
    \item 
    If the entries of \code{x} are not all integers, then
    \code{rounding(x)} returns the smallest number of digits \code{d}
    after the decimal point
    such that \code{\link[base]{round}(x, digits=d)} is identical to
    \code{x}.
    For example if \code{rounding(x) = 2} then the entries of
    \code{x} are rounded to 2 decimal places, and are multiples of 0.01.
    \item 
    If all the entries of \code{x} are integers, then
    \code{rounding(x)} returns \code{-d}, where
    \code{d} is the smallest number of digits \emph{before} the decimal point
    such that \code{\link[base]{round}(x, digits=-d)} is identical to
    \code{x}.
    For example if \code{rounding(x) = -3} then the entries of
    \code{x} are multiples of 1000.
    If \code{rounding(x) = 0} then the entries of \code{x} are integers
    but not multiples of 10.
    \item
    If all entries of \code{x} are equal to 0, a value of 0 is returned.
  }

}
\value{
  An integer.
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{round.ppp} in package \code{spatstat.geom}.
}
\examples{
   rounding(c(0.1, 0.3, 1.2))
   rounding(c(1940, 1880, 2010))
   rounding(0)
}
\keyword{math}
