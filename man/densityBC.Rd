\name{densityBC}
\alias{densityBC}
\title{Kernel Density Estimation with Optional Boundary Correction}
\description{
  Fixed-bandwidth kernel density
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
  \item{x}{Numeric vector of observed values.}
  \item{kernel}{String specifying kernel.
    Options are
    \code{"gaussian"}, \code{"rectangular"},
    \code{"triangular"},
    \code{"epanechnikov"},
    \code{"biweight"},
    \code{"cosine"} and \code{"optcosine"}.
    (Partial matching is used).
    Options are described in the help for \code{\link[stats]{density.default}}.
  }
  \item{bw,h}{
    Alternative specifications of the scale factor for the kernel.
    The bandwidth \code{bw} is the standard deviation of the
    kernel (this agrees with the argument \code{bw} in
    \code{\link[stats]{density.default}}).
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
    \code{"none"}, \code{"weighted"}, \code{"convolution"},
    \code{"reflection"}, \code{"bdrykern"} and \code{"JonesFoster"}. 
  }
  \item{fast}{
    Logical value specifying whether to perform the calculation rapidly
    using the Fast Fourier Transform (\code{fast=TRUE})
    or to use slower, exact code (\code{fast=FALSE}, the default).
  }
  \item{internal}{
    Internal use only.
  }
  \item{\dots}{
    Additional arguments are ignored.
  }
}
\details{
  This function computes a fixed-bandwidth kernel estimate of
  a probability density on the real line, or the positive half-line,
  including optional boundary corrections
  for truncation of the density onto the positive half line.

  Weighted estimates are supported, including negative weights.
  Weights are not renormalised to sum to 1. The resulting
  probability density estimate is not renormalised to integrate to 1.

  Options for the smoothing kernel
  are described in the help for \code{\link[stats]{density.default}}.
  The default is the Epanechnikov (truncated quadratic) kernel.
  
  If \code{zerocor} is missing, or given as \code{"none"},
  this function computes the fixed-bandwidth kernel estimator of the
  probability density on the real line,
  using \code{\link[stats]{density.default}}.
  The estimated probability density (unnormalised) is
  \deqn{
    \widehat f(x) = \sum_{i=1}^n w_i \; \kappa(x - x_i)
  }{
    f(x) = sum[i=1,...,n] w[i] * kappa(x - x[i])
  }
  where \eqn{x_1,\ldots,x_n}{x[1], ..., x[n]} are the data values,
  \eqn{w_1,\ldots,w_n}{w[1], ..., w[n]} are the weights (defaulting
  to \eqn{w_i = 1/n}{w[i] = 1/n}),
  and \eqn{\kappa}{kappa} is the kernel, a probability density on the
  real line. 

  If \code{zerocor} is given, the probability density is assumed to be
  confined to the positive half-line; the numerical values in \code{x}
  must all be non-negative; and a boundary correction is
  applied to compensate for bias arising due to truncation at the origin:
  \describe{
    \item{\code{zerocor="weighted"}:}{
      The contribution from each data point \eqn{x_i}{x[i]}
      is weighted by the factor \eqn{1/m(x_i)}{1/m(x[i])}
      where \eqn{m(x) = 1 - F(-x)} is the total mass of the kernel centred on
      \eqn{x} that lies in the positive half-line, and \eqn{F(x)} is the
      cumulative distribution function of the kernel.
      The corrected estimate is
      \deqn{
	\widehat f_W(x) = \sum_{i=1}^n w_i \; \frac{\kappa(x - x_i)}{1-F(-x_i)}
      }{
	fW(x) = sum[i=1,...,n] w[i] * kappa(x - x[i])/(1-F(-x[i]))
      }
      This is the \dQuote{cut-and-normalization} method of
      Gasser and \ifelse{latex}{\out{M\"{u}ller}}{Mueller} (1979).
      Effectively the kernel is renormalized so that it integrates to 1,
      and the adjusted kernel conserves mass.
    }
    \item{\code{zerocor="convolution"}:}{
      The estimate of the density \eqn{f(x)} is
      weighted by the factor \eqn{1/m(x)} where \eqn{m(r) = 1 - F(-x)}
      is given above.
      The corrected estimate is
      \deqn{
	\widehat f_C(x) = \sum_{i=1}^n w_i \; \frac{\kappa(x - x_i)}{1-F(-x)}
      }{
	fC(x) = sum[i=1,...,n] w[i] * kappa(x - x[i])/(1-F(-x))
      }
      This is the \dQuote{convolution}, \dQuote{uniform}
      or \dQuote{zero-order} boundary correction method
      often attributed to Diggle (1985).
      This correction does not conserve mass.
      It is faster to compute than the weighted correction.
    }
    \item{\code{zerocor="reflection"}:}{
      if the kernel centred at data point \eqn{x_i}{x[i]}
      has a tail that lies on the negative half-line, this tail is
      reflected onto the positive half-line.
      The corrected estimate is
      \deqn{
	\widehat f_R(x) = \sum_{i=1}^n w_i \; 
	              [ \kappa(x - x_i) + \kappa(-x - x_i) ]
      }{
	fR(x) = sum[i=1,...,n] w[i] * ( kappa(x - x[i]) + kappa(-x - x[i]) )
      }
      This is the \dQuote{reflection} method first proposed by
      Boneva et al (1971). This correction conserves mass.
      The estimated density always has zero derivative at the origin,
      \eqn{\widehat f_R^\prime(0) = 0}{fR'(0) = 0}, which may or
      may not be desirable.
    }
    \item{\code{zerocor="bdrykern"}:}{
      The density estimate is computed using the
      Linear Boundary Kernel associated with the chosen kernel
      (Wand and Jones, 1995, page 47).
      The estimated (unnormalised) probability density is
      \deqn{
	\widehat f_B(x) =  \sum_{i=1}^n w_i \;
	[ A(x) + (x-x_i) B(x)] \kappa(x - x_i)
      }{
	fB(x) = sum[i=1,...,n] w[i] * 
	( A(x) + (x-x[i]) B(x)) * kappa(x - x[i])
      }
      where \eqn{A(x) = a_2(x)/D(x)}{A(x) = a[2](x)/D(x)} and
      \eqn{B(x) = -a_1(x)/D(x)}{B(x) = -a[1](x)/D(x)}
      with
      \eqn{D(x) = a_0(x) a_2(x) - a_1(x)^2}{D(x) = a[0](x) a[2](x) - a[1](x)^2}
      where
      \eqn{
	a_k(x) = \int_{-\infty}^x t^k \kappa(t) dt.
      }{
	a[k](x) = integral[-Inf,x] ( t^k * kappa(t) dt).
      }
      That is, when estimating the density \eqn{f(x)} for values of
      \eqn{x} close to zero (defined as \eqn{x < h} for all kernels
      except the Gaussian), the kernel contribution
      \eqn{k_h(x - x_i)}{k[h](x - x[i])} is multiplied by a
      term that is a linear function of \eqn{x - x_i}{x - x[i]},
      with coefficients depending on \eqn{x}.
      This correction does not conserve mass and may result in
      negative values, but is asymptotically optimal.
      Computation time for this estimate is greater than for
      the options above.
    }
    \item{\code{zerocor="JonesFoster"}:}{
      The modification of the Boundary Kernel estimate
      proposed by Jones and Foster (1996) is computed:
      \deqn{
	\widehat f_{JF}(x) =
	\widehat f_C(x)
	\exp\left( \frac{\widehat f_B(x)}{\widehat f_C(r)} - 1 \right)
      }{
	fJF(x) = fC(x) exp(fB(x)/fC(x) - 1)
      }
      where \eqn{\widehat f_C(r)}{fC(r)} is the convolution estimator
      and \eqn{\widehat f_B(r)}{fB(r)} is the linear boundary kernel estimator.
      This ensures that the estimate is always nonnegative
      and retains the asymptotic optimality of the linear boundary
      kernel.
      Computation time for this estimate
      is greater than for all the options above.
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
\seealso{
  \code{\link[stats]{density.default}}.
  
  \code{\link{dkernel}} for the kernel itself.

  \code{\link{densityAdaptiveKernel.default}} for adaptive
  (variable-bandwidth) estimation.
}
\references{
  \ournewpaper

  Boneva, L.I., Kendall, D.G. and Stefanov, I. (1971)
  Spline transformations: three new diagnostic aids for the
  statistical data-analyst (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{33}, 1--70.

  Diggle, P.J. (1985)
  A kernel method for smoothing point process data.
  \emph{Journal of the Royal Statistical Society, Series C (Applied Statistics)}, 
  \bold{34} 138--147.
  
  Gasser, Th. and \ifelse{latex}{\out{M\"{u}ller}}{Mueller}, H.-G. (1979).
  Kernel estimation of regression functions.
  In Th. Gasser and M. Rosenblatt (editors)
  \emph{Smoothing Techniques for Curve Estimation}, pages
  23--68. Springer, Berlin.

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

