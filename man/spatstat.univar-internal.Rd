\name{spatstat.univar-internal}
\title{Internal Functions}
\alias{as.breakpts}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{check.hist.lengths}
\alias{km.rs.opt}
\alias{make.even.breaks}
\alias{match.kernel}
\alias{print.ewcdf}
\alias{quantilefun.interpolatedCDF}
\alias{IdenticalRowPair}
\alias{IdenticalRows}
\alias{StieltjesCalc}
\alias{StieltjesCalc.stepfun}
\description{
  Undocumented Functions
}
\usage{
as.breakpts(\dots)
breakpts(val, maxi, even, npos, step)
breakpts.from.r(r)
check.hist.lengths(hist, breaks)
make.even.breaks(bmax, npos, bstep)
km.rs.opt(o, cc, d, breaks, KM, RS)
make.even.breaks(bmax, npos, bstep)
match.kernel(kernel)
\method{print}{ewcdf}(x, digits, \dots)
\method{quantilefun}{interpolatedCDF}(x, \dots, type)
IdenticalRowPair(i,j,a,b)
IdenticalRows(i,j,a,b)
StieltjesCalc(M, f, \dots)
\method{StieltjesCalc}{stepfun}(M, f, \dots)
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
