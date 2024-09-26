\name{bw.taylor}
\alias{bw.taylor}
\title{
  Bandwidth Selection for Kernel Density Estimation by Non-Random Bootstrap
}
\description{
  Use Taylor's non-random bootstrap technique to select
  the bandwidth for kernel density estimation on the real line.
}
\usage{
bw.taylor(x, \dots, srange = NULL, useC = TRUE)
}
\arguments{
  \item{x}{
    Numeric vector.
  }
  \item{\dots}{
    Ignored.
  }
  \item{srange}{
    Range of bandwidths to be considered.
    A numeric vector of length 2.
  }
  \item{useC}{
    Logical value specifying whether to use faster C code.
  }
}
\details{
  This function selects a bandwidth for kernel density estimation
  of a probability density on the real line,
  using the numeric data \code{x} and assuming a Gaussian kernel. The result
  is the numeric value of the standard deviation of the Gaussian kernel.

  The function uses the method of Taylor (1989) who showed that,
  when using the Gaussian kernel, the optimisation criterion
  can be computed rapidly from the data without any randomised
  resampling. 

  The domain of the probability density is assumed to be
  the entire real line. Boundary correction is not currently implemented.
  % (but modifications could
  % theoretically be made to allow this in exchange for an increase in
  % computational cost).

  The result of \code{bw.taylor} is a single numeric value
  giving the selected bandwidth.
%  If \code{ns} is missing or \code{NULL}, then the optimal value
%  is found rapidly using \code{\link[stats]{optimise}}.
%  If \code{ns} is an integer, then the optimisation criterion is
%  evaluated at \code{ns} candidate values of bandwidth, and the result
%  also belongs to the class \code{"bw.optim"}; it contains the
%  values of the optimisation criterion at all candidate bandwidths,
%  and can be plotted as shown in the examples.
}
\value{
  A single numeric value.
%  If \code{ns} is given, the result also belongs to the class
%  \code{"bw.optim"}, and contains the values of the optimisation
%  criterion, which can be plotted using \code{plot}.
}
\references{
  Taylor, C.C. (1989) Choice of the Smoothing Parameter in Kernel
  Density Estimation, \emph{Biometrika} \bold{76} 4, 705--712. 
}
\author{
  \tilman and \adrian.
}
\seealso{
  \code{\link[stats]{bw.nrd}} in the \pkg{stats} package
  for standard bandwidth selectors.
}
\examples{
  x <- rnorm(30)
  bw.taylor(x)
}
\keyword{univar}
\keyword{smooth}
\concept{Bandwidth selection}
