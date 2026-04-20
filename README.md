# spatstat.univar

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.univar)](http://CRAN.R-project.org/package=spatstat.univar) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.univar)](https://github.com/spatstat/spatstat.univar)

This member of the `spatstat` package family
contains code for estimation of **univariate probability distributions**.

You are viewing the GitHub repository which holds
the latest **development version** of `spatstat.univar`.
For the latest public release on CRAN, click the green badge above.

### Overview 

`spatstat.univar` contains code for
estimation of univariate probability distributions, including:

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
- *kernel density estimation*
for one dimensional probability densities, including
fixed- and variable-bandwidth kernel estimators,
and boundary corrections for densities restricted to the positive half-line;
- *kernels* 
including calculation of the probability density, cumulative distribution
function, quantiles, random generation, moments and partial
moments of the standard smoothing kernels;
- *special distributions:*
including calculation of the one-dimensional heat kernel in an interval,
and the distribution of a weighted sum of chi-squared variables;
- *integration:*
Numerical integration including Stieltjes integrals
and indefinite integrals.

Some code for discrete probability distributions is provided by 
[`spatstat.random`](https://github.com/spatstat/spatstat.random).

### Installing the package

This repository contains the _development version_ of
`spatstat.univar`. The easiest way to install the development version
is to start R and type

```R
repo <- c('https://spatstat.r-universe.dev', 'https://cloud.r-project.org')
install.packages("spatstat.univar", dependencies=TRUE, repos=repo)
```

To install the latest _public release_ of `spatstat.univar`,
type

```R
install.packages("spatstat.univar")
```

