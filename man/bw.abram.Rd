\name{bw.abram}
\alias{bw.abram}
\title{
  Abramson's Adaptive Bandwidths
}
\description{
  Computes adaptive smoothing bandwidths
  according to the inverse-square-root rule of Abramson (1982).
}
\usage{
bw.abram(X, h0, \dots)
}
\arguments{
  \item{X}{
    Data to be smoothed.
  }
  \item{h0}{
    Global smoothing bandwidth. A numeric value.
  }
  \item{\dots}{
    Additional arguments passed to methods.
  }
}
\details{
  This function computes adaptive smoothing bandwidths for a dataset,
  using the methods of Abramson (1982) and Hall and Marron (1988).

  The function \code{bw.abram} is generic.
  There is a default method \code{\link{bw.abram.default}}.
  The \pkg{spatstat} package family includes methods for spatial objects.
}
\seealso{
  \code{\link{bw.abram.default}}
}
\value{
  See the documentation for the particular method.
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
  \adrian.
}
\keyword{nonparametric}
\concept{Adaptive smoothing}
\concept{Bandwidth selection}
