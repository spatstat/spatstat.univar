\name{densityAdaptiveKernel}
\alias{densityAdaptiveKernel}
\title{
  Adaptive Kernel Estimation of Density or Intensity
}
\description{
  Computes an adaptive estimate of probability density or intensity
  using a variable-bandwidth smoothing kernel.
}
\usage{
  densityAdaptiveKernel(X, \dots)
}
\arguments{
  \item{X}{
    Data to be smoothed. 
  }
  \item{\dots}{
    Additional arguments passed to methods.
  }
}
\details{
  This generic function computes an adaptive kernel estimate of 
  probability density or intensity.

  The function \code{densityAdaptiveKernel} is generic.
  The \pkg{spatstat} package family includes methods for spatial objects.
}
\value{
  See documentation for each method.
}
\author{
  \adrian and \tilman.
}
\seealso{
  \code{\link{bw.abram}}.
}
\keyword{nonparametric}
