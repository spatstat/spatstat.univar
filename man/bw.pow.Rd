\name{bw.pow}
\alias{bw.pow}
\title{
  Variable Bandwidths Proportional to a Power of the Data Value
}
\description{
  Computes variable smoothing bandwidths intended to be proportional to
  the observed data values, raised to a given power.
}
\usage{
 bw.pow(X, h0, POW = 0.75, trim = 5, \dots)
}
\arguments{
  \item{X}{
    Data for which bandwidths should be calculated. A numeric vector
    of positive values.
  }
  \item{h0}{
    A scalar value giving the global smoothing bandwidth
    in the same units as \code{X}. The default is
    \code{\link[stats]{bw.nrd0}(X)}.
  }
  \item{POW}{
   Numeric value. The exponent of the power transformation to
   be applied to \code{X}. 
  }
  \item{trim}{
    A trimming value required to curb excessively large bandwidths.
    See Details. The default is sensible in most cases.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This function computes adaptive smoothing bandwidths
  for the data values in \code{X}.
  Larger data values are assigned larger bandwidths.

  Bandwidths are proportional to \code{X^POW}.
  The bandwidth at location \eqn{u} is
  \deqn{
    h(u) = \mbox{\texttt{h0}}
    * \mbox{min}[ \frac{u^{\mbox{\texttt{POW}}}}{\gamma}, \mbox{\texttt{trim}} ]
  }{
    h(u) = h0 * min(u^POW/\gamma, trim)
  }
  where \eqn{\gamma} is the geometric mean of the values
  \eqn{u^{\mbox{\texttt{POW}}}}{u^POW}.
  This allows the global bandwidth \code{h0} to be considered on
  the same scale as a corresponding fixed bandwidth.
}
\value{
 A numeric vector of the same length as \code{X}.
}
\seealso{
  \code{\link[spatstat.univar]{bw.abram}},
  \code{\link[stats]{bw.nrd0}}.
}
% \references{
% Brainpower.
% }
\author{
  \tilman. Adapted by \adrian.
}
\examples{
  xx <- sort(rexp(10))
  bb <- bw.pow(xx)
  signif(rbind(xx, bb), 3)
}
\keyword{nonparametric}

