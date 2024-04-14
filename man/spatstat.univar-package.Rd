\name{spatstat.univar-package}
\alias{spatstat.univar-package} 
\alias{spatstat.univar} 
\docType{package}
\title{The spatstat.univar Package}
\description{
  The \pkg{spatstat.univar} package
  belongs to the \pkg{spatstat} family of packages.
  It provides utilities for estimating the probability distribution
  of one-dimensional (real-valued) data.
}
\details{
  This package is a member of the \pkg{spatstat}
  family of packages. It provides utilities for
  estimation of the probability distribution of
  one-dimensional (i.e. numerical, real-valued) data.
  The utilities include:
  \describe{
    \item{kernel density estimation:}{
      including
      unnormalised weighted densities, and cumulative distribution
      functions of density estimates.
    }
    \item{weighted distributions and weighted statistics:}{
      including weighted empirical cumulative distributions, weighted median,
      weighted quantiles, calculating the CDF from a density estimate
    }
    \item{estimation for right-censored data:}{
      including Kaplan-Meier, reduced-sample and other estimators
      of the cumulative distribution function and hazard function
      from right-censored data
    }
    \item{quantiles:}{
      including calculation of quantiles from an empirical cumulative
      distribution or a kernel density estimate
    }
    \item{kernels:}{
      including calculation of the probability density, cumulative distribution
      function, quantiles, random generation, moments and partial
      moments of the standard smoothing kernels
    }
    \item{heat kernel:}{
      calculation of the one-dimensional heat kernel in an interval
    }
    \item{integration:}{
      Numerical integration including Stieltjes integrals
      and indefinite integrals.
    }
  }

  The facilities are described in more detail below.
  
  \bold{Kernel density estimation}

  The package supports 
  fixed-bandwidth and variable-bandwidth kernel estimation
  of probability densities from numerical data.
  It provides boundary corrections for kernel estimates
  of densities on the positive half-line (applicable when the original
  observations are positive numbers) for both fixed-bandwidth
  and variable-bandwidth estimates.

  If the observations have numerical weights associated with them,
  these weights will not be automatically normalised, and indeed
  the weights may be negative or zero. This is 
  unlike the standard \R method \code{\link[stats]{density.default}}.

  The main functions are:
  \tabular{ll}{
    \code{\link{unnormdensity}}
    \tab
    extension of \code{\link[stats]{density.default}}
    allowing weights to be negative or zero.
    \cr
    \code{\link{CDF.density}} \tab cumulative distribution function from
    kernel density estimate
  }

  \bold{Weighted distributions and weighted statistics}

  Weighted versions of standard operations such as the histogram
  and empirical distribution function are provided:
  
  \tabular{ll}{
    \code{\link{whist}} \tab weighted histogram \cr
    \code{\link{ewcdf}} \tab weighted empirical cumulative distribution
    function \cr
    \code{\link{mean.ewcdf}} \tab mean of weighted ecdf \cr
    \code{\link{quantile.ewcdf}} \tab quantiles of weighted ecdf \cr
    \code{\link{weighted.median}} \tab weighted median of numeric values \cr
    \code{\link{weighted.quantile}} \tab weighted quantile of numeric values \cr
  }

  \bold{Estimation for right-censored data}

  Facilities are provided for estimating the probability distribution
  of right-censored lifetimes (non-negative real random variables).
  
  \tabular{ll}{
    \code{\link{kaplan.meier}} \tab Kaplan-Meier estimator of cumulative
    distribution function and hazard rate, from right-censored data \cr
    \code{\link{reduced.sample}} \tab reduced-sample estimator of
    cumulative distribution function, from right-censored data
  }
    
  \bold{Quantiles}

  Facilities are provided for computing the quantiles
  of a probability distribution, given estimates of the
  probability density or the cumulative distribution function
  and so on.
  
  \tabular{ll}{
    \code{\link{CDF.density}} \tab cumulative distribution function from
    kernel density estimate \cr
    \code{\link{quantile.density}} \tab quantiles of kernel density estimate \cr
    \code{\link{quantile.ewcdf}} \tab quantiles of weighted ecdf \cr
    \code{\link{quantilefun}} \tab quantiles as a function \cr
    \code{\link{quantilefun.ewcdf}} \tab quantiles as a function \cr
    \code{\link{weighted.quantile}} \tab weighted quantile of numeric values \cr
    \code{\link{transformquantiles}} \tab transform the quantiles of a dataset
  }

  \bold{Kernels}

  The standard \R function \code{\link[stats]{density.default}}
  recognises a list of smoothing kernels by name:
  \code{"gaussian"}, \code{"rectangular"}, \code{"triangular"},
  \code{"epanechnikov"}, \code{"biweight"}, \code{"cosine"}
  and \code{"optcosine"}. For these kernels, \pkg{spatstat.univar}
  provides various characteristics:
  \tabular{ll}{
    \code{\link{dkernel}}
    \tab
    probability density of the kernel
    \cr
    \code{\link{pkernel}}
    \tab
    cumulative distribution function of the kernel
    \cr
    \code{\link{qkernel}}
    \tab
    quantiles of the kernel
    \cr
    \code{\link{rkernel}}
    \tab
    generate simulated realisations from the kernel
    \cr
    \code{\link{kernel.factor}}
    \tab
    scale factor relating bandwidth to half-width of kernel
    \cr
    \code{\link{kernel.moment}}
    \tab
    partial moment of kernel
    \cr
    \code{\link{kernel.squint}}
    \tab
    integral of squared kernel
  }

  \bold{Heat kernels}

  The heat kernel in an interval can be calculated.

  \tabular{ll}{  
    \code{\link{hotrod}}
    \tab
    calculate the heat kernel in an interval
  }

  \bold{Integration}

  A few facilities are provided for calculating integrals
  of real functions.
  
  \tabular{ll}{
    \code{\link{indefinteg}} \tab indefinite integral \cr
    \code{\link{integral.density}} \tab integral of a kernel density
    estimate \cr
    \code{\link{stieltjes}} \tab Stieltjes integral
  }

}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \adrian, \tilman, \martinH, \ege, \rolf and Greg McSwiggan.
}
\keyword{package}
\keyword{distribution}
\keyword{nonparametric}
\keyword{smooth}

