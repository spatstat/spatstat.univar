\name{dsocs}
\alias{dsocs}
\alias{psocs}
\alias{qsocs}
\alias{rsocs}
\title{
  Distribution of a Linear Combination of Chi-Squared Variables
}
\description{
  Probability density, cumulative probability distribution function,
  quantile function, and random generation for a linear combination of
  \eqn{\chi^2}{chi^2} random variables.
}
\usage{
dsocs(x, lambda, log= FALSE, method=c("Wood", "Farebrother"))
psocs(q, lambda, lower.tail=TRUE, log.p=FALSE, method=c("Wood", "Farebrother"))
qsocs(p, lambda, lower.tail=TRUE, log.p=FALSE, method=c("Wood", "Farebrother")) 
rsocs(n, lambda, approx=FALSE)
}
\arguments{
  \item{x,q}{Numeric vector of quantiles}
  \item{p}{Numeric vector of probabilities}
  \item{n}{Integer. Number of random values to generate.}
  \item{lambda}{Numeric vector of coefficients of the linear combination.}
  \item{lower.tail}{
    Logical value. If \code{TRUE} (the default), probabilities
    are \eqn{P[X \le x]}{P[X <= x]}; otherwise, \eqn{P[X > x]}/
  }
  \item{log,log.p}{
    Logical value. If \code{TRUE}, probabilities are given as natural
    logarithms \eqn{\log(p)}{log(p)}.
  }
  \item{approx}{
    Logical value specifying whether to generate random values
    from the approximate distribution (\code{approx=TRUE})
    or according to the exact distribution (\code{approx=FALSE}, the
    default).
  }
  \item{method}{
    Character string (partially matched) specifying
    whether to use Wood's approximation (\code{method="Wood"}, the
    default) or Farebrother's algorithm (\code{method="Farebrother"}).
  }
}
\details{
  These functions calculate the probability density (\code{dsocs}),
  cumulative distribution function (\code{psocs}) and
  quantile function (\code{qsocs}), and generate random realisations
  (\code{rsocs}) of, a random variable \eqn{Y} such that
  \deqn{
    Y = \sum_{j=1^m} \lambda_j X_j^2
  }{
    Y = sum[j] lambda[j] * X[j]^2
  }
  where \eqn{X_1,\ldots,X_m}{X[1],...,X[m]} are independent, identically
  distributed standard normal random variables, and
  \eqn{\lambda_1,\ldots,\lambda_m}{lambda[1], ..., lambda[m]} are
  nonnegative coefficients. Thus, \eqn{Y} is a linear combination of
  independent \eqn{\chi^2}{chi^2} random variables
  each having 1 degree of freedom.

  The argument \code{lambda} contains the coefficients
  \eqn{\lambda_1,\ldots,\lambda_m}{lambda[1], ..., lambda[m]}.

  \itemize{
    \item if \code{method="Wood"} (the default),
    the approximation method of Wood (1989) is used.
    The distribution of \eqn{Y} is approximated by an
    \eqn{F} distribution, except in marginal cases when it is approximated
    by a gamma distribution or inverse-gamma distribution.
    \item if \code{method="Farebrother"}, the
    numerical algorithm of Farebrother (1984) is used.
    The cumulative distribution function and the probability density
    are evaluated to a high level of precision
    using a truncated infinite series. The quantile function is
    evaluated by root-finding.
  }

  For the random generator \code{rsocs}, the default behaviour is to
  generate realisations of \eqn{Y} according to its exact distribution
  by constructing \eqn{\chi^2}{chi^2} random variables and summing.
  If \code{approx=TRUE}, then the realisations are generated from the
  (Wood's method) approximate distribution instead.
}
\value{
  A numeric vector.
}
\references{
  Farebrother, R.W. (1984)
  Algorithm AS 204:
  The distribution of a positive linear combination of \eqn{\chi^2}{chi^2}
  random variables. \emph{Applied Statistics} \bold{33} #3 (1984) 332--339.

  Wood, A.T.A (1989)
  An F approximation to the distribution of a linear combination
  of chi-squared variables.
  \emph{Communications in Statistics -- Simulation and Computation}
  \bold{18}:4, 1439--1456.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats:FDist]{pf}},
  \code{\link[stats:GammaDist]{pgamma}}
}
\examples{
   N <- if(interactive()) 1000 else 50
   co <- c(1, 0.5, 0.3)
   xx <- rsocs(N, co)
   plot(density(xx))
   curve(dsocs(x, co), add=TRUE, col=4)
   qqplot(xx, qsocs(ppoints(N), co))
   abline(0,1)
   plot(ppoints(N), psocs(sort(xx), co))
   abline(0,1)
}
\keyword{distribution}
\keyword{univar}
\keyword{datagen}
