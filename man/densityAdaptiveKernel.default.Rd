\name{densityAdaptiveKernel.default}
\alias{densityAdaptiveKernel.default}
\title{
  Adaptive Kernel Estimation of Probability Density
}
\description{
  Computes an adaptive estimate of probability density
  from numeric data, using a variable-bandwidth smoothing kernel.
}
\usage{
\method{densityAdaptiveKernel}{default}(X, bw, \dots,
       weights = NULL,
       zerocor=c("none", "weighted", "convolution",
                 "reflection", "bdrykern", "JonesFoster"),
       at = c("grid", "data"), ngroups=Inf, fast=TRUE)
}
\arguments{
  \item{X}{
    Data to be smoothed. A numeric vector.
  }
  \item{bw}{
    Smoothing bandwidths. Either a numeric vector of the same length as
    \code{X} giving the bandwidth associated with each data value,
    or a \code{function} in the \R language that provides the smoothing
    bandwidth at any desired location. The default is to compute
    bandwidths using \code{\link[spatstat.univar]{bw.abram.default}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[stats]{density.default}}
    controlling the range of \code{x} values at which the density must
    be estimated, when \code{at="grid"}.
  }
  \item{weights}{
    Optional. Numeric vector of weights attached to each value in \code{X}.
  }
  \item{zerocor}{
    Character string (partially matched) specifying a boundary
    correction. This is appropriate when \code{X} contains only positive
    values.
  }
  \item{at}{
    String (partially matched) specifying whether to evaluate the
    probability density only at the data points (\code{at="data"})
    or on a grid of \code{x} values (\code{at="grid"}, the default).
  }
  \item{ngroups}{
    Integer, \code{Inf} or \code{NULL}.
    If \code{ngroups = Inf}, the density estimate will be computed
    exactly using C code.
    If \code{ngroups} is finite, then the fast subdivision technique
    of Davies and Baddeley (2018) will be applied.
    If \code{ngroups = NULL} then a default rule is used to choose
    an efficient number of groups.
  }
  \item{fast}{
    Logical value specifying whether to use the Fast Fourier Transform
    to accelerate computations, when appropriate.
  }
}
\details{
  This function computes an adaptive kernel estimate of 
  probability density on the real line (if \code{zerocor="none"})
  or on the positive real line (if \code{zerocor} is another value).

  The argument \code{bw} specifies the smoothing bandwidths to be
  applied to each of the points in \code{X}. It may be a numeric vector
  of bandwidth values, or a function yielding the
  bandwidth values.

  If the values in \code{X} are \eqn{x_1,\ldots,x_n}{x[1], ..., x[n]}
  and the corresponding bandwidths are
  \eqn{\sigma_1,\ldots,\sigma_n}{\sigma[1], ..., \sigma[n]}
  then the adaptive kernel estimate of intensity at a location \eqn{u} is
  \deqn{
    \hat\lambda(u) = \sum_{i=1}^n k(u, x_i, \sigma_i)
  }{
    \lambda(u) = sum[i] e(x[i], k(u, x[i], \sigma[i])
  }
  where \eqn{k(u, v, \sigma)} is the value at \eqn{u}
  of the (possibly edge-corrected) smoothing kernel with bandwidth \eqn{\sigma}
  induced by a data point at \eqn{v}.

  Exact computation of the estimate above can be time-consuming:
  it takes \eqn{n} times longer than fixed-bandwidth smoothing.
  
  The partitioning method of Davies and Baddeley (2018)
  accelerates this computation by partitioning the range of
  bandwidths into \code{ngroups} intervals,
  correspondingly subdividing \code{X} into
  \code{ngroups} subsets according to bandwidth,
  and applying fixed-bandwidth smoothing to each subset.

  If \code{ngroups=NULL} then we use a default rule where \code{ngroups}
  is the integer part of the square root of
  the number of points in \code{X}, so that the computation time is
  only about \eqn{\sqrt{n}}{sqrt(n)} times slower than fixed-bandwidth
  smoothing. Any positive value of \code{ngroups}
  can be specified by the user. Specifying \code{ngroups=Inf} enforces exact
  computation of the estimate without partitioning. Specifying
  \code{ngroups=1} is the same as fixed-bandwidth smoothing with
  bandwidth \code{sigma=median(bw)}.
}
\section{Bandwidths and Bandwidth Selection}{
  The function \code{densityAdaptiveKernel.default}
  computes one adaptive estimate of probability density,
  determined by the smoothing bandwidth values \code{bw}.

  Typically the bandwidth values are computed by first computing
  a pilot estimate of the intensity, then using \code{\link[spatstat.univar]{bw.abram.default}}
  to compute the vector of bandwidths according to Abramson's rule.
  This involves specifying a global bandwidth \code{h0}.
}
\value{
  If \code{at="data"}, a numeric vector of the same length as \code{X}.
  If \code{at="grid"}, a probability density object of class \code{"density"}.
}
\references{
  Davies, T.M. and Baddeley, A. (2018)
  Fast computation of spatially adaptive kernel estimates.
  \emph{Statistics and Computing}, \bold{28}(4), 937-956.\cr

  Hall, P. and Marron, J.S. (1988)
  Variable window width kernel density estimates of probability
  densities.
  \emph{Probability Theory and Related Fields}, \bold{80}, 37-49.\cr

  Silverman, B.W. (1986)
  \emph{Density Estimation for Statistics and Data Analysis}.
  Chapman and Hall, New York.
}
\author{
  \adrian and \tilman.
}
\seealso{
  \code{\link[spatstat.univar]{bw.abram.default}}
}
\examples{
  xx <- rexp(100, rate=5)
  plot(density(xx))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
  plot(densityAdaptiveKernel(xx, at="grid"))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
  plot(densityAdaptiveKernel(xx, at="grid", zerocor="w"))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
  plot(densityAdaptiveKernel(xx, at="grid", zerocor="c"))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
  plot(densityAdaptiveKernel(xx, at="grid", zerocor="r"))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
  plot(densityAdaptiveKernel(xx, at="grid", zerocor="b"))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
  plot(densityAdaptiveKernel(xx, at="grid", zerocor="J"))
  curve(5 * exp(-5 * x), add=TRUE, col=3)
}
\keyword{nonparametric}
