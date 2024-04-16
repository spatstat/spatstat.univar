\name{spatstat.univar-internal}
\title{Internal Functions}
\alias{as.breakpts}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{check.hist.lengths}
\alias{make.even.breaks}
\alias{match.kernel}
\alias{print.ewcdf}
\alias{StieltjesCalc}
\alias{StieltjesCalc.stepfun}
\alias{quantilefun.interpolatedCDF}
\description{
  Undocumented Functions
}
\usage{
as.breakpts(\dots)
breakpts(val, maxi, even, npos, step)
breakpts.from.r(r)
check.hist.lengths(hist, breaks)
make.even.breaks(bmax, npos, bstep)
match.kernel(kernel)
\method{print}{ewcdf}(x, digits, \dots)
StieltjesCalc(M, f, \dots)
\method{StieltjesCalc}{stepfun}(M, f, \dots)
\method{quantilefun}{interpolatedCDF}(x, \dots, type)
}
\details{
  These functions are not documented.
}
\value{
  \code{as.breakpts}, \code{breakpts}, \code{breakpts.from.r} and
  \code{make.even.breaks} return an object of class \code{"breakpts"}
  which defines a sequence of breakpoints, identifies whether they are
  equally-spaced and so on. \code{match.kernel} returns a single
  character string. \code{check.hist.lengths} returns \code{NULL}.
  \code{quantilefun.interpolatedCDF} returns a function.
  The \code{print} method returns \code{NULL}. Other functions are
  not documented.
}
\keyword{internal}
