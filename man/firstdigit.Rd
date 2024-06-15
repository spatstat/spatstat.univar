\name{firstdigit}
\alias{firstdigit}
\alias{lastdigit}
\alias{ndigits}
\title{
  Digits in Decimal Representation
}
\description{
  Find the first or last digit in the decimal representation of a number.
}
\usage{
firstdigit(x)
lastdigit(x)
ndigits(x)
}
\arguments{
  \item{x}{A numeric value or numeric vector.}
}
\details{
  \code{firstdigit(x)} finds the first (most significant) digit,
  \code{lastdigit(x)} finds the last (least significant) digit,
  and \code{ndigits(x)} finds the number of digits,
  in the decimal representation of each entry of \code{x}.
  The decimal representation is truncated at the number of digits
  available for double precision numbers on the hardware, usually 15.
}
\value{
  An integer or integer vector of the same length as \code{x}.
}
\author{
  \adrian
}
\seealso{
  \code{\link[spatstat.univar]{rounding}}
}
\examples{
firstdigit(42)
lastdigit(42)
ndigits(42)

firstdigit(-0.1234)
lastdigit(-0.1234)
ndigits(-0.1234)

firstdigit(0)
lastdigit(0)
ndigits(0)
}
\keyword{math}