\name{densityBC}
\alias{densityBC}
\title{Kernel Density Estimation with Optional Boundary Correction}
\description{
  A simple implementation of fixed-bandwidth kernel density
  estimation on the real line, or the positive real half-line,
  including optional corrections for a boundary at zero. 
}
\usage{
  densityBC(x, kernel = "epanechnikov", bw=NULL,
      \dots,
      h=NULL,
      adjust = 1,
      weights = rep(1, length(x))/length(x), from, to = max(x), n = 256,
      zerocor = c("none", "weighted", "convolution", "reflection",
                  "bdrykern", "JonesFoster"),
      fast=FALSE,
      internal=list())
}
\arguments{
  \item{x}{Numeric vector.}
  \item{kernel}{String specifying kernel.
    Options are
    \code{"gaussian"}, \code{"rectangular"},
    \code{"triangular"},
    \code{"epanechnikov"},
    \code{"biweight"},
    \code{"cosine"} and \code{"optcosine"}.
    (Partial matching is used).
  }
  \item{bw,h}{
    Alternative specifications of the scale factor for the kernel.
    The bandwidth \code{bw} is the standard deviation of the
    kernel (this agrees with the argument \code{bw} in
    \code{\link[stats]{density.default}}.
    The rescale factor \code{h} is the factor by which
    the `standard form' of the kernel is rescaled. For the Epanechnikov
    kernel, \code{h = bw * sqrt(5)} is the
    half-width of the support, while for the Gaussian kernel,
    \code{h = bw} is the standard deviation.
    Either \code{bw} or \code{h} should be given,
    and should be a single numeric value,
    or a character string indicating a bandwidth selection rule
    as described in \code{\link[stats]{density.default}}.
  }
  \item{adjust}{
    Numeric value used to rescale the bandwidth \code{bw} and halfwidth
    \code{h}. The bandwidth used is \code{adjust * bw}. This makes it
    easy to specify values like \sQuote{half the default} bandwidth.
  }
  \item{weights}{
    Numeric vector of weights associated with \code{x}.
    The weights are not required to sum to 1, and will not be normalised
    to sum to 1. The weights may include negative values.
  }
  \item{from,to}{
    Lower and upper limits of interval on which density should be
    computed.
    The default value of \code{from} is \code{from=min(x)} if
    \code{zerocor="none"}, and \code{from=0} otherwise.
  }
  \item{n}{
    Number of \eqn{r} values for which density should be computed.
  }
  \item{zerocor}{
    String (partially matched) specifying a correction for the boundary effect
    bias at \eqn{r=0} when estimating a density on the positive
    half line. Possible values are
    \code{"none"}, \code{"weighted"}, \code{"convolution"}, \code{"reflection"}
    and \code{"bdrykern"}. 
  }
  \item{fast}{
    Logical value specifying whether to perform the calculation rapidly
    using the Fast Fourier Transform (\code{fast=TRUE})
    or to use slower, exact code (\code{fast=FALSE}, the default).
    Option \code{zerocor="bdrykern"} is not available when \code{fast=TRUE}.
  }
  \item{internal}{
    Internal use only.
  }
  \item{\dots}{
    Additional arguments are ignored.
  }
}
\details{
  If \code{zerocor} is absent or given as \code{"none"},
  this function computes the fixed bandwidth kernel estimator of the
  probability density on the real line. 

  If \code{zerocor} is given, it is assumed that the density
  is confined to the positive half-line, and a boundary correction is
  applied:
  \describe{
    \item{weighted}{The contribution from each point \eqn{x_i}{x[i]}
      is weighted by the factor \eqn{1/m(x_i)}{1/m(x[i])}
      where \eqn{m(x) = 1 - F(-x)} is the total mass of the kernel centred on
      \eqn{x} that lies in the positive half-line, and \eqn{F(x)} is the
      cumulative distribution function of the kernel}
    \item{convolution}{The estimate of the density \eqn{f(r)} is
      weighted by the factor \eqn{1/m(r)} where \eqn{m(r) = 1 - F(-r)}
      is given above.
    }
    \item{reflection}{
      if the kernel centred at data point \eqn{x_i}{x[i]}
      has a tail that lies on the negative half-line, this tail is
      reflected onto the positive half-line.
    }
    \item{bdrykern}{The density estimate is computed using the
      Boundary Kernel associated with the chosen kernel
      (Wand and Jones, 1995, page 47).
      That is, when estimating the density \eqn{f(r)} for values of
      \eqn{r} close to zero (defined as \eqn{r < h} for all kernels
      except the Gaussian), the kernel contribution
      \eqn{k_h(r - x_i)}{k[h](r - x[i])} is multiplied by a
      term that is a linear function of \eqn{r - x_i}{r - x[i]}.
    }
    \item{JonesFoster}{
      The modification of the Boundary Kernel estimate
      proposed by Jones and Foster (1996), equal to
      \eqn{
	\overline f(r) \exp( \hat f(r)/\overline f(r) - 1)
      }{
	f#(r) exp(f*(r)/f#(r) - 1)
      }
      where \eqn{\overline f(r)}{f#(r)} is the convolution estimator
      and \eqn{\hat f(r)}{f*(r)} is the boundary kernel estimator.
    }
  }
  If \code{fast=TRUE}, the calculations are performed rapidly using 
  \code{\link[stats]{density.default}} which employs the Fast Fourier
  Transform. If \code{fast=FALSE} (the default), the calculations are
  performed exactly using slower C code. 
}
\value{
  An object of class \code{"density"} as described in the help file
  for \code{\link[stats]{density.default}}. It contains at least the entries
  \item{x }{Vector of \eqn{x} values}
  \item{y }{Vector of density values \eqn{y= f(x)}}
}
\examples{
  sim.dat <- rexp(500)
  fhatN <- densityBC(sim.dat, "biweight", h=0.4)
  fhatB <- densityBC(sim.dat, "biweight", h=0.4, zerocor="bdrykern")
  plot(fhatN, ylim=c(0,1.1), main="density estimates")
  lines(fhatB, col=2)
  curve(dexp(x), add=TRUE, from=0, col=3)
  legend(2, 0.8,
     legend=c("fixed bandwidth", "boundary kernel", "true density"),
     col=1:3, lty=rep(1,3))

}
\references{
  \ournewpaper

  Jones, M.C. and Foster, P.J. (1996)
  A simple nonnegative boundary correction method for kernel density
  estimation.
  \emph{Statistica Sinica}, \bold{6} (4) 1005--1013.

  Wand, M.P. and Jones, M.C. (1995)
  \emph{Kernel Smoothing}.
  Chapman and Hall.
}
\author{\adrian and \martinH.}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}

