\name{kernel.moment}
\alias{kernel.moment}
\title{Incomplete Moment of Smoothing Kernel}
\description{
  Computes the complete or incomplete \eqn{m}th moment of a
  smoothing kernel.
}
\usage{
  kernel.moment(m, r, kernel = "gaussian", mean=0, sd=1/kernel.factor(kernel))
}
\arguments{
  \item{m}{
    Exponent (order of moment).
    An integer.
  }
  \item{r}{
    Upper limit of integration for the incomplete moment.
    A numeric value or numeric vector.
    Set \code{r=Inf} to obtain the complete moment.
  }
  \item{kernel}{
    String name of the kernel.
    Options are
    \code{"gaussian"}, \code{"rectangular"},
    \code{"triangular"},
    \code{"epanechnikov"},
    \code{"biweight"},
    \code{"cosine"} and \code{"optcosine"}.
    (Partial matching is used).
  }
  \item{mean,sd}{
    Optional numerical values giving the mean and standard deviation
    of the kernel.
  }
}
\details{
  Kernel estimation of a probability density in one dimension
  is performed by \code{\link[stats]{density.default}}
  using a kernel function selected from the list above.
  For more information about these kernels,
  see \code{\link[stats]{density.default}}. 

  The function \code{kernel.moment} computes the integral 
  \deqn{
    \int_{-\infty}^r t^m k(t) dt
  }{
    integral[-Inf][r] t^m k(t) dt
  }
  where \eqn{k(t)} is the selected kernel, \eqn{r} is the upper limit of
  integration, and \eqn{m} is the exponent or order.

  Note that, if \code{mean} and \code{sd} are not specified, the
  calculations assume that \eqn{k(t)} is the \bold{standard form} of the kernel,
  which has support \eqn{[-1,1]} and
  standard deviation \eqn{sigma = 1/c} where \code{c = kernel.factor(kernel)}.

  The code uses the explicit analytic expressions when
  \code{m = 0, 1, 2} and numerical integration otherwise.
}
\value{
  A single number, or a numeric vector of the same length as \code{r}.
}
\seealso{
  \code{\link[stats]{density.default}},
  \code{\link{dkernel}},
  \code{\link{kernel.factor}},
  \code{\link{kernel.squint}}
}
\examples{
   kernel.moment(1, 0.1, "epa")
   curve(kernel.moment(2, x, "epa"), from=-1, to=1)
}
\author{
  \adrian
  and \martinH.
}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}
