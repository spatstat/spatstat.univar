\name{bw.abram.default}
\alias{bw.abram.default}
\title{
  Abramson's Adaptive Bandwidths For Numeric Data
}
\description{
  Computes adaptive smoothing bandwidths for numeric data,
  according to the inverse-square-root rule of Abramson (1982).
}
\usage{
 \method{bw.abram}{default}(X, h0, \dots,
    at = c("data", "grid"),
    pilot = NULL, hp = h0, trim = 5, smoother = density.default)
}
\arguments{
  \item{X}{
    Data for which bandwidths should be calculated. A numeric vector.
  }
  \item{h0}{
    A scalar value giving the global smoothing bandwidth
    in the same units as \code{X}.
    The default is \code{h0=\link[stats]{bw.nrd0}(X)}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[stats]{density.default}}
    (or to \code{smoother}) controlling the range of values \code{x}
    at which the density will be estimated, when \code{at="grid"}.
  }
  \item{at}{
    Character string (partially matched) specifying whether to
    compute bandwidth values only at the data points of \code{X}
    (\code{at = 'data'}, the default) or on a grid of \code{x} values
    (\code{at = 'grid'}).
  }
  \item{pilot}{
    Optional. Specification of a pilot density (possibly unnormalised).
    Either a numeric vector giving the pilot density for each
    data point of \code{X}, a \code{function} in the \R language,
    or a probability density estimate (object of class \code{"density"}).
    If \code{pilot=NULL} the pilot density
    is computed by applying fixed-bandwidth density estimation to
    \code{X} using bandwidth \code{hp}.
  }
  \item{hp}{
    Optional. A scalar pilot bandwidth, used for estimation
    of the pilot density, if \code{pilot} is not given.
  }
  \item{trim}{
    A trimming value required to curb excessively large bandwidths.
    See Details. The default is sensible in most cases.
  }
  \item{smoother}{
    Smoother for the pilot.
    A function or character string, specifying the function
    to be used to compute the pilot estimate when
    \code{pilot} is \code{NULL}.
  }
}
\details{
  This function computes adaptive smoothing bandwidths
  using the methods of Abramson (1982) and Hall and Marron (1988).

  The function \code{\link[spatstat.univar]{bw.abram}} is generic. The function
  \code{bw.abram.default} documented here is the default method
  which is designed for numeric data.
  
  If \code{at="data"} (the default) a smoothing bandwidth is
  computed for each data point in \code{X}. Alternatively if
  \code{at="grid"} a smoothing bandwidth is computed for
  a grid of \code{x} values.

  Under the Abramson-Hall-Marron rule, the bandwidth at location \eqn{u} is
  \deqn{
    h(u) = \mbox{\texttt{h0}}
    * \mbox{min}[ \frac{\tilde{f}(u)^{-1/2}}{\gamma}, \mbox{\texttt{trim}} ]
  }{
    h(u) = h0 * min(\tilde{f}(u)^{-1/2}/\gamma, trim)
  }
  where \eqn{\tilde{f}(u)} is a pilot estimate of the 
  probability density. The variable bandwidths are rescaled by \eqn{\gamma}, the
  geometric mean of the \eqn{\tilde{f}(u)^{-1/2}} terms evaluated at the
  data; this allows the global bandwidth \code{h0} to be considered on
  the same scale as a corresponding fixed bandwidth. The trimming value
  \code{trim} has the same interpretation as the required `clipping' of
  the pilot density at some small nominal value (see Hall and Marron,
  1988), to necessarily prevent extreme bandwidths (which
  can occur at very isolated observations).

  The pilot density or intensity is determined as follows:
  \itemize{
    \item If \code{pilot} is a \code{function} in the \R language,
    this is taken as the pilot density.
    \item If \code{pilot} is a probability density estimate
    (object of class \code{"density"} produced by
    \code{\link[stats]{density.default}}) then this is taken as the 
    pilot density.
    \item If \code{pilot} is \code{NULL}, then the pilot intensity is
    computed as a fixed-bandwidth kernel
    intensity estimate using \code{\link{density.default}} applied to
    the data \code{X} using
    the pilot bandwidth \code{hp}.
  }
  In each case the pilot density is renormalised to
  become a probability density, and then the Abramson rule is applied.

  Instead of calculating the pilot as a fixed-bandwidth density
  estimate, the user can specify another density estimation procedure
  using the argument \code{smoother}. This should be either a function
  or the character string name of a function. It will replace
  \code{\link{density.default}} as the function used to calculate the
  pilot estimate. The pilot estimate will be computed as
  \code{smoother(X, sigma=hp, ...)} if \code{pilot} is \code{NULL},
  or \code{smoother(pilot, sigma=hp, ...)} if \code{pilot} is a point
  pattern. If \code{smoother} does not recognise the argument name
  \code{sigma} for the smoothing bandwidth, then \code{hp} is effectively
  ignored.
}
\value{
  Either a numeric vector of the same length as \code{X}
  giving the Abramson bandwidth for each point
  (when \code{at = "data"}, the default),
  or a \code{function} giving the Abramson bandwidths
  as a function of location.
}
\seealso{
  \code{\link[spatstat.univar]{bw.abram}},
  \code{\link[stats]{bw.nrd0}}.
}
\references{
  Abramson, I. (1982)
  On bandwidth variation in kernel estimates --- a square root law.
  \emph{Annals of Statistics}, \bold{10}(4), 1217-1223. 

  Hall, P. and Marron, J.S. (1988)
  Variable window width kernel density estimates of probability
  densities.
  \emph{Probability Theory and Related Fields}, \bold{80}, 37-49. 

  Silverman, B.W. (1986)
  \emph{Density Estimation for Statistics and Data Analysis}.
  Chapman and Hall, New York.
}
\author{
  \tilman. Adapted by \adrian.
}
\examples{
  xx <- rexp(20)
  bw.abram(xx)
}
\keyword{nonparametric}

