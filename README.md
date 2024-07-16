# spatstat.univar

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.univar)](http://CRAN.R-project.org/package=spatstat.univar) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.univar)](https://github.com/spatstat/spatstat.univar)

This is the newest member of the `spatstat` package family.
It contains code for estimation of univariate probability distributions,
including:

- *weighted distributions and weighted statistics*
including weighted empirical cumulative distributions, weighted median,
weighted quantiles, calculating the CDF from a density estimate;
- *estimation for right-censored data* 
including Kaplan-Meier, reduced-sample and other estimators
of the cumulative distribution function and hazard function
from right-censored data;
- *quantiles* 
including calculation of quantiles from an empirical cumulative
distribution or a kernel density estimate;
- *kernels* 
including calculation of the probability density, cumulative distribution
function, quantiles, random generation, moments and partial
moments of the standard smoothing kernels;
- *heat kernel:*
calculation of the one-dimensional heat kernel in an interval;
- *integration:*
Numerical integration including Stieltjes integrals
and indefinite integrals.

Some of the code has been extracted from `spatstat.geom`,
`spatstat.random` and `spatstat.explore`, while some is new.
