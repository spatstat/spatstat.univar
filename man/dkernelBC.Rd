\name{dkernelBC}
\alias{dkernelBC}
\title{
  Boundary-corrected Kernel Density Function
}
\description{
  Computes the boundary-corrected version of a smoothing kernel density
  function.
}
\usage{
dkernelBC(x, mean, sd = 1, kernel = "gaussian",
       zerocor = c("none", "weighted", "convolution",
                    "reflection", "bdrykern"))
}
\arguments{
  \item{x}{
    Numeric. Values of the function argument, at which the
    function should be evaluated.
  }
  \item{mean}{
    Numeric. The mean of the uncorrected kernel.
  }
  \item{sd}{
    Numeric value. The standard deviation of the uncorrected kernel.
  }
  \item{kernel}{
    Character string giving the name of the kernel
    as recognised by \code{\link[spatstat.univar]{match.kernel}}.
}
  \item{zerocor}{
    String (partially matched) specifying a correction for the boundary effect
    bias at \eqn{r=0} when estimating a density on the positive
    half line. Possible values are
    \code{"none"}, \code{"weighted"}, \code{"convolution"}, \code{"reflection"},
    and \code{"bdrykern"}.
  }
}
\details{
  The kernel density function identified by \code{kernel}
  with standard deviation \code{sd} and mean \code{mean}
  will be computed, and truncated onto the positive half-line.
  The boundary correction specified by \code{zerocor} will then
  be applied. The result is the vector of corrected density values.
}
\value{
  Numeric value or numeric vector.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{densityBC}} to compute a density estimate using
  the boundary-corrected kernel.

  \code{\link[spatstat.univar]{dkernel}} to compute the un-corrected
  kernel density function, and \code{\link[stats]{density.default}} to compute
  an uncorrected density estimate.

  \code{\link[spatstat.univar]{match.kernel}} for the list of
  recognised names of kernels.
}
\examples{
  curve(dkernelBC(x, mean=1, zerocor="none"), to=5)
  curve(dkernelBC(x, mean=1, zerocor="weighted"), to=5)
  curve(dkernelBC(x, mean=1, zerocor="reflection"), to=5)
  curve(dkernelBC(x, mean=1, zerocor="convolution"), to=5)
  curve(dkernelBC(x, mean=1, zerocor="bdrykern"), to=5)
}
