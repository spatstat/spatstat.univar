\name{weighted.median}
\alias{weighted.median}
\alias{weighted.quantile}
\alias{weighted.var}
\title{
  Weighted Median, Quantiles or Variance
}
\description{
  Compute the median, quantiles or variance of a set of numbers which have
  weights associated with them.
}
\usage{
weighted.median(x, w, na.rm = TRUE, type=2, collapse=TRUE)

weighted.quantile(x, w, probs=seq(0,1,0.25), na.rm = TRUE, type=4, collapse=TRUE)

weighted.var(x, w, na.rm = TRUE)
}
\arguments{
  \item{x}{
    Data values.
    A vector of numeric values, for which the median or quantiles are required.
  }
  \item{w}{
    Weights.
    A vector of nonnegative numbers, of the same length as \code{x}.
  }
  \item{probs}{
    Probabilities for which the quantiles should be computed.
    A numeric vector of values between 0 and 1.
  }
  \item{na.rm}{
    Logical. Whether to ignore \code{NA} values.
  }
  \item{type}{
    Integer specifying the rule for calculating the median or quantile,
    corresponding to the rules available for
    \code{\link[stats]{quantile}}.
    The currently available choices are \code{type=1}, \code{2} or \code{4}.
    See Details.
  }
  \item{collapse}{Logical value specifying whether
    duplicated values in \code{x} should be pooled (replacing them by
    a unique \code{x} value whose weight is the sum of the associated
    weights).
  }
}
\details{
  The \code{i}th observation \code{x[i]} is treated as having
  a weight proportional to \code{w[i]}.

  The weighted median is a value \code{m}
  such that the total weight of data less than or equal to \code{m}
  is equal to half the total weight. More generally, the weighted quantile with
  probability \code{p} is a value \code{q} 
  such that the total weight of data less than or equal to \code{q}
  is equal to \code{p} times the total weight.

  If one of the data values does satisfy this requirement, then
  \itemize{
    \item if \code{type=1} or \code{type=4}, this value is returned;
    \item if \code{type=2}, 
    the average of this value and the next largest value is returned.
  }

  If there is no such value, then 
  \itemize{
    \item if \code{type=1} or \code{type=2}, the next largest value is returned;
    \item if \code{type=4}, the result is a value obtained by
    linearly interpolating between the next smallest and next largest values.
  }

  Type 1 is the right-continuous inverse of the weighted
  cumulative distribution function. 

  Type 2 (with \code{collapse=FALSE}) is consistent with the
  traditional definition of the sample median.

  Note that the default settings are different for \code{weighted.median}
  and \code{weighted.quantile}.
}
\value{
  \code{weighted.median} returns a numeric value.
  \code{weighted.quantile} returns a numeric vector
  of the same length as \code{probs}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats]{quantile}}, \code{\link[stats]{median}}.
}
\examples{
  x <- 1:20
  w <- runif(20)
  weighted.median(x, w)
  weighted.quantile(x, w)
  weighted.var(x, w)
}
\keyword{math}
