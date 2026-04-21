\name{farebro}
\alias{farebro}
\title{
  Distribution of a Weighted Sum of Noncentral Chi-Squared Variables
}
\description{
  Calculates the cumulative probability distribution function
  and the probability density function of 
  a weighted sum of noncentral chi-squared variables.
}
\usage{
farebro(lambda, mult = 1, delta = 0, x, \dots,
        maxit = 1e+06, eps = 1e-06, warn = TRUE)
}
\arguments{
  \item{lambda}{
    Numeric vector of positive weights.
  }
  \item{mult}{
    Integer value or integer vector of multiplicities (degrees of freedom).
  }
  \item{delta}{
    Numeric value or numeric vector of noncentrality parameters.
  }
  \item{x}{
    Numeric value or vector giving the function argument
    (quantiles for which the cumulative probability and density are desired).
  }
  \item{\dots}{
    Ignored.
  }
  \item{maxit}{
    Integer. Maximum number of iterations.
  }
  \item{eps}{
    Positive numeric value. Maximum error in the probability value. 
  }
  \item{warn}{
    Logical value specifying whether to issue a warning
    if an error code is returned.
  }
}
\details{
  This function implements the algorithm of
  Farebrother (1984) for computing the cumulative distribution
  function \eqn{F(x) = P(X \le x)}{F(x) = P(X <= x)}
  and probability density function \eqn{f(x) = F^\prime(x)}{f(x) = F'(x)}
  of the random variable
  \deqn{
    Y = \sum_{j=1}^n \lambda_j X_j
  }{
    Y = sum[j=1,..,n] lambda[j] * X[j]
  }
  where \eqn{X_1,\ldots,X_n}{X[1], ..., X[n]} are independent
  random variables and \eqn{X_j}{X[j]} has a noncentral
  \eqn{\chi^2}{chi^2} distribution with \eqn{m_j}{m[j]} degrees of
  freedom and noncentrality parameter \eqn{\delta_j}{delta[j]}.

  The arguments \code{lambda}, \code{mult} and \code{delta}
  should be vectors of equal length, or will be recycled so that
  they have equal length.

  The algorithm uses an exact mathematical representation of \eqn{F(x)}
  as an infinite series, and a mathematical upper bound on the error
  of the truncated series, to compute an approximation to \eqn{F(x)}
  which is guaranteed to be accurate to within the specified error
  \code{eps}. The algorithm involves calculating the probability density
  \eqn{f(x)} as well. It returns the values of both \eqn{F(x)} and
  \eqn{f(x)}.

  An error code \code{ifault} is returned for each entry of \code{x}.
  Nonzero error codes are translated to text warnings if
  \code{warn=TRUE}.
}
\value{
  A data frame, with one row for each entry in \code{x},
  and columns \code{x} (input value),
  \code{p} (cumulative probability),
  \code{d} (probability density)
  and \code{ifault} (error code). 
}
\references{
  Farebrother, R.W. (1984)
  Algorithm AS 204:
  The distribution of a positive linear combination of \eqn{\chi^2}{chi^2}
  random variables. \emph{Applied Statistics} \bold{33} (3) 332--339.
}
\author{
  R.W. Farebrother, adapted by \adrian.
  Original \proglang{Algol} code of Farebrother (1984)
  translated to \proglang{C} by Adrian Baddeley
  and modified to handle a vector of quantile values.
  \R interface by Adrian Baddeley.
}
\seealso{
  \code{\link{dsocs}}
}
\examples{
  ## Farebrother (1984) Table 1, example Q_5 (rows 13-15)
  farebro(c(7,3), c(6,2), c(6,2), c(20, 100, 200))
}
\keyword{distribution}
\keyword{univar}
