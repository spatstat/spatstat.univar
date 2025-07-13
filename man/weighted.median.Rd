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
weighted.median(x, w, na.rm = TRUE, type=2, collapse=FALSE)

weighted.quantile(x, w, probs=seq(0,1,0.25), na.rm = TRUE, type=4, collapse=FALSE)

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
    If \code{w} is missing or \code{NULL}, the default weights
    are all equal to 1.
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
    The currently available choices are \code{type=1}, \code{2},
    \code{3} and \code{4}.
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

  The weighted sample median is a value \code{m} 
  such that the total weight of data less than or equal to \code{m}
  is equal to half the total weight. More generally,
  the weighted sample quantile with
  probability \code{p} is a value \code{q} 
  such that the total weight of data less than or equal to \code{q}
  is equal to \code{p} times the total weight.

  Define the weighted empirical cumulative distribution function
  \deqn{
    F(x) = \sum_{i: x_i \le x} w_i / \sum_{i=1}^n w_i
  }{
    F(x) = sum (w_i | x_i <= x) / sum(w_i)
  }
  That is, \eqn{F(x)} is the fraction of total weight that is
  associated with data values \eqn{x_i}{x[i]}
  less than or equal to the value \eqn{x}.

  The weighted quantile for probability \eqn{p} is a number \eqn{q}
  defined so that \eqn{F(q) = p} wherever possible.
  There are different definitions of the quantile depending on
  how this should be achieved.
  
  For unweighted data, 
  there are 9 different definitions of the sample median and sample
  quantile, which enjoy slightly different properties. These 9 different
  definitions are explained in the help for
  \code{\link[stats]{quantile.default}}. The user's choice of algorithm is
  selected using the argument \code{type}.
  
  For weighted data, the first
4
  of the 9 definitions of sample quantile
  have been generalised to weighted quantiles.
  The functions \code{weighted.median} and \code{weighted.quantile}
  documented here provide these definitions of weighted sample
  quantile. The user's choice of algorithm is again
  selected using the argument \code{type}.
  
  Suppose the data values have been arranged in increasing order
  as
  \eqn{x_{[1]} \le x_{[2]} \le \ldots \le x_{[n]}}{x[1] <= x[2] <= ... x[n]}.
  If one of the data values \eqn{x_{[k]}}{x[k]} satisfies
  \eqn{F(x_{[k]})=p}{F(x[k]) = p} exactly, then
  \itemize{
    \item if \code{type=1}, \code{type=3} or \code{type=4},
    the code returns this value, \eqn{x_{[k]}}{x[k]};
    \item if \code{type=2},
    the code returns
    the average of this value and the next largest value,
    \eqn{(x_{[k]} + x_{[k+1]})/2}.
  }

  If there is no data value satisfying
  \eqn{F(x_{[k]})=p}{F(x[k]) = p} exactly, then
  the code finds the two adjacent values
  \eqn{x_{[k]}}{x[k]} and \eqn{x_{[k+1]}}{x[k+1]} which satisfy
  \eqn{F(x_{[k]}) < p}{F(x[k]) < p} and
  \eqn{F(x_{[k+1]}) > p}{F(x[k+1]) > p}, and defines the quantile
  as follows:
  
  \itemize{
    \item if \code{type=1} or \code{type=2}, the code returns
    the larger value \eqn{x_{[k+1]}}{x[k+1]};
    \item if \code{type=3}, the code returns the value
    which minimises the discrepancy, that is,
    if we define
    \eqn{d_k = |F(x_{[k]}) - p|}{d[k] = |F(x[k]) - p|}
    and 
    \eqn{d_{k+1} = |F(x_{[k+1]}) - p|}{d[k+1] = |F(x[k+1]) - p|},
    then 
    \itemize{
      \item
      if \eqn{d_k < d_{k+1}}{d[k] < d[k+1]}, 
      the code returns \eqn{x_{[k]}}{x[k]};
      \item
      if \eqn{d_k > d_{k+1}}{d[k] > d[k+1]}, 
      the code returns \eqn{x_{[k+1]}}{x[k+1]};
      \item
      if \eqn{d_k = d_{k+1}}{d[k] = d[k+1]}, 
      then the even-numbered value
      is returned, that is, the code returns \eqn{x_{[k^\ast]}}{x[k*]} where
      \eqn{k^\ast}{k*} equals \eqn{k} if \eqn{k} is even, and equals
      \eqn{k+1} if \eqn{k+1} is even.
    }
    \item if \code{type=4}, the code returns
    a value obtained by linearly interpolating between the
    two adjacent values \eqn{x_{[k]}}{x[k]} and \eqn{x_{[k+1]}}{x[k+1]},
    that is, it returns the value
    \eqn{(1-a) x_{[k]} + a x_{[k+1]}}{(1-a) * x[k] + a * x[k+1]}
    where
    \eqn{a = (p - F(x_{[k]}))/(F(x_{[k+1]}) - F(x_{[k]}))}{a = (p -
      F(x[k]))/(F(x[k+1]) - F(x[k])}.
  }

  For very small probabilities \eqn{p < F(x_{[1]})}{p < F(x[1])} the value
  \eqn{x_{[1]}}{x[1]} is returned.
  For very large probabilities \eqn{p > F(x_{[n]})}{p > F(x[n])} the value
  \eqn{x_{[n]}}{x[n]} is returned.  

  Type 1 is the left-continuous quantile function.

  Type 2 is consistent with the
  traditional definition of the sample median.

  Types 1 and 3 always return a value selected from the input data \code{x},
  while types 2 and 4 sometimes return values that are interpolated
  between the input data values.

  Note that the default settings are different for \code{weighted.median}
  and \code{weighted.quantile}.

%  The implementations of type 3 and 5 are
  The implementation of type 3 is
  experimental and may be changed.
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
  weighted.quantile(x, w, probs=(0:4)/4)
  weighted.var(x, w)
}
\keyword{math}
