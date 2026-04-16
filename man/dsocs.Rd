\name{dsocs}
\alias{dsocs}
\alias{psocs}
\alias{qsocs}
\alias{rsocs}
\title{
  Approximate Distribution of a Linear Combination of Chi-Squared Variables
}
\description{
  Approximations to the probability density, cumulative
  probability distribution function, quantile function,
  and random generation for a linear combination of
  \eqn{\chi^2}{chi^2} random variables.
}
\usage{
dsocs(x, lambda, log= FALSE)
psocs(q, lambda, lower.tail=TRUE, log.p=FALSE)
qsocs(p, lambda, lower.tail=TRUE, log.p=FALSE) 
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
}
\details{
  These functions calculate approximations to
  the probability density (\code{dsocs}),
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
  
  The approximation method of Wood (1989) is used.
  The distribution of \eqn{Y} is approximated by an
  \eqn{F} distribution, except in marginal cases when it is approximated
  by a gamma distribution or inverse-gamma distribution.

  For the random generator \code{rsocs}, the default behaviour is to
  generate realisations of \eqn{Y} according to its exact distribution
  by constructing \eqn{\chi^2}{chi^2} random variables and summing.
  If \code{approx=TRUE}, then the realisations are generated from the
  approximate distribution instead.
}
\value{
  A numeric vector.
}
\references{
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
